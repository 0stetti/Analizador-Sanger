# ==========================================
# REQUISITOS (guardar num ficheiro requirements.txt):
# streamlit>=1.30.0
# biopython>=1.81
# plotly>=5.18.0
# ==========================================

import streamlit as st
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import plotly.graph_objects as go
import io

# Configuração inicial da página
st.set_page_config(page_title="Sanger Pro Viewer", layout="wide")

# --- FUNÇÕES DE PROCESSAMENTO ---
def obter_dados_sanger(record):
    """
    Extrai os traços (ondas), a sequência original e a posição exata dos picos (PLOC).
    """
    # Canais padrão ABI para G, A, T, C
    canais = ['DATA9', 'DATA10', 'DATA11', 'DATA12'] 
    mapa_bases = {'DATA9': 'G', 'DATA10': 'A', 'DATA11': 'T', 'DATA12': 'C'}
    cores = {'G': 'black', 'A': 'green', 'T': 'red', 'C': 'blue'}
    
    tracos = {}
    
    # Extrair os dados brutos (raw data)
    if 'abif_raw' in record.annotations:
        raw = record.annotations['abif_raw']
        
        # 1. Extrair as ondas de fluorescência
        for canal in canais:
            if canal in raw:
                base = mapa_bases[canal]
                tracos[base] = list(raw[canal])
        
        # 2. Extrair as Posições dos Picos (Peak Locations - PLOC)
        # O PLOC diz-nos exatamente em que coordenada X a máquina detetou a base
        ploc = raw.get('PLOC_1', raw.get('PLOC_2', []))
        
        return tracos, ploc, cores
    return None, None, None

# --- INTERFACE PRINCIPAL ---
st.title("🧬 Sanger Pro: Visualização e Alinhamento")

# --- BARRA LATERAL (CONTROLES E UPLOAD) ---
with st.sidebar:
    st.header("1. Carregar Dados")
    ficheiro_carregado = st.file_uploader("Ficheiro .ab1", type=["ab1"])
    
    st.markdown("---")
    st.header("2. Controlos de Visualização")
    escala_vertical = st.slider("Amplitude dos Picos (Zoom Vertical)", 1.0, 10.0, 1.0, 0.1)
    
    st.markdown("---")
    st.header("3. Referência (Opcional)")
    seq_referencia = st.text_area("Sequência Teórica (ex: mclover 3)", height=150).upper().strip()

# --- ÁREA DE TRABALHO ---
if ficheiro_carregado:
    try:
        # Ler o ficheiro binário em memória
        dados_bytes = ficheiro_carregado.read()
        record = SeqIO.read(io.BytesIO(dados_bytes), "abi")
        
        # Extrair os dados complexos do cromatograma
        tracos, plocs, cores = obter_dados_sanger(record)
        
        # Inicializar o estado da sessão para manter as edições do utilizador
        if 'seq_editada' not in st.session_state or st.session_state.get('id_ficheiro') != ficheiro_carregado.name:
            st.session_state['seq_editada'] = str(record.seq)
            st.session_state['id_ficheiro'] = ficheiro_carregado.name

        st.success(f"Ficheiro '{ficheiro_carregado.name}' carregado com sucesso!")

        # Separar a interface em separadores (Tabs)
        tab_grafico, tab_alinhamento = st.tabs(["📊 Cromatograma e Edição", "🔍 Alinhamento"])

        # ==========================================
        # SEPARADOR 1: GRÁFICO E EDIÇÃO
        # ==========================================
        with tab_grafico:
            st.subheader("Cromatograma Interativo")
            
            fig = go.Figure()
            valor_maximo = 0

            # 1. Desenhar as ondas
            if tracos:
                for base, dados in tracos.items():
                    # Aplicar a escala vertical para aumentar picos baixos
                    dados_escalados = [d * escala_vertical for d in dados]
                    if dados_escalados:
                        valor_maximo = max(valor_maximo, max(dados_escalados))
                    
                    fig.add_trace(go.Scatter(
                        y=dados_escalados,
                        name=base,
                        mode='lines',
                        line=dict(color=cores[base], width=1),
                        hoverinfo='skip' # Melhora a performance
                    ))

            # 2. Desenhar as letras no topo dos picos (usando o PLOC)
            seq_atual = st.session_state['seq_editada']
            limite = min(len(plocs), len(seq_atual))
            
            fig.add_trace(go.Scatter(
                x=list(plocs)[:limite], 
                y=[valor_maximo * 1.05] * limite, # Colocar a 5% acima do pico mais alto
                text=list(seq_atual)[:limite],
                mode="text",
                textfont=dict(size=14, color="black"),
                name="Bases Chamadas"
            ))

            # 3. Configurar o layout com a barra de deslocamento (Range Slider)
            fig.update_layout(
                height=450,
                showlegend=True,
                plot_bgcolor='white',
                margin=dict(l=10, r=10, t=30, b=10),
                xaxis=dict(
                    title="Posição do Traço",
                    rangeslider=dict(visible=True), # Barra de deslocamento horizontal
                    showgrid=True,
                    gridcolor='lightgrey'
                ),
                yaxis=dict(
                    title="Intensidade",
                    showgrid=True,
                    gridcolor='lightgrey',
                    fixedrange=False # Permite zoom no eixo Y
                )
            )

            st.plotly_chart(fig, use_container_width=True)

            # 4. Caixa de Edição
            st.markdown("### ✏️ Editor de Bases")
            st.info("Altera a sequência abaixo. O gráfico será atualizado automaticamente com as novas letras.")
            
            col1, col2 = st.columns([4, 1])
            with col1:
                nova_seq = st.text_area(
                    "Sequência Extraída",
                    value=st.session_state['seq_editada'],
                    height=150
                ).upper().strip()
                
                # Se o utilizador editar, atualizar o estado e recarregar
                if nova_seq != st.session_state['seq_editada']:
                    st.session_state['seq_editada'] = nova_seq
                    st.rerun()
            
            with col2:
                st.write("**Estatísticas**")
                st.metric("Tamanho Original", len(plocs))
                st.metric("Tamanho Editado", len(nova_seq))
                if len(nova_seq) != len(plocs):
                    st.warning("⚠️ O tamanho foi alterado. O alinhamento visual com os picos pode perder a sincronia.")

        # ==========================================
        # SEPARADOR 2: ALINHAMENTO
        # ==========================================
        with tab_alinhamento:
            st.subheader("Alinhamento com Sequência Teórica")
            
            if seq_referencia:
                if st.button("Executar Alinhamento", type="primary"):
                    # Configurar o algoritmo de alinhamento
                    aligner = PairwiseAligner()
                    aligner.mode = 'local'
                    aligner.match_score = 2
                    aligner.mismatch_score = -1
                    aligner.open_gap_score = -2
                    aligner.extend_gap_score = -1
                    
                    # Executar usando a sequência EDITADA pelo utilizador
                    alinhamentos = aligner.align(seq_referencia, st.session_state['seq_editada'])
                    melhor_alinhamento = alinhamentos[0]
                    
                    st.metric("Pontuação do Alinhamento (Score)", melhor_alinhamento.score)
                    st.text("Visão do Alinhamento:")
                    st.code(str(melhor_alinhamento), language='text')
            else:
                st.warning("Insere uma sequência de referência na barra lateral para efetuar o alinhamento.")

    except Exception as e:
        st.error(f"Ocorreu um erro ao processar o ficheiro: {e}")

else:
    # Ecrã inicial quando não há ficheiros
    st.info("👈 Começa por carregar um ficheiro .ab1 na barra lateral.")
