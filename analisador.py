# ==========================================
# REQUISITOS (guardar num ficheiro requirements.txt):
# streamlit>=1.35.0   <-- ATENÇÃO: É necessário o Streamlit 1.35 ou superior para eventos de clique!
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
    canais = ['DATA9', 'DATA10', 'DATA11', 'DATA12'] 
    mapa_bases = {'DATA9': 'G', 'DATA10': 'A', 'DATA11': 'T', 'DATA12': 'C'}
    cores = {'G': 'black', 'A': 'green', 'T': 'red', 'C': 'blue'}
    
    tracos = {}
    
    if 'abif_raw' in record.annotations:
        raw = record.annotations['abif_raw']
        
        for canal in canais:
            if canal in raw:
                base = mapa_bases[canal]
                tracos[base] = list(raw[canal])
        
        ploc = raw.get('PLOC_1', raw.get('PLOC_2', []))
        return tracos, ploc, cores
    return None, None, None

# --- INTERFACE PRINCIPAL ---
st.title("🧬 Sanger Pro: Visualização e Alinhamento")

with st.sidebar:
    st.header("1. Carregar Dados")
    ficheiro_carregado = st.file_uploader("Ficheiro .ab1", type=["ab1"])
    
    st.markdown("---")
    st.header("2. Controlos de Visualização")
    escala_vertical = st.slider("Amplitude dos Picos (Zoom Vertical)", 1.0, 10.0, 1.0, 0.1)
    
    st.markdown("---")
    st.header("3. Referência (Opcional)")
    seq_referencia = st.text_area("Sequência Teórica (ex: mclover 3)", height=150).upper().strip()

if ficheiro_carregado:
    try:
        dados_bytes = ficheiro_carregado.read()
        record = SeqIO.read(io.BytesIO(dados_bytes), "abi")
        
        tracos, plocs, cores = obter_dados_sanger(record)
        
        if 'seq_editada' not in st.session_state or st.session_state.get('id_ficheiro') != ficheiro_carregado.name:
            st.session_state['seq_editada'] = str(record.seq)
            st.session_state['id_ficheiro'] = ficheiro_carregado.name

        tab_grafico, tab_alinhamento = st.tabs(["📊 Cromatograma Interativo", "🔍 Alinhamento Global"])

        with tab_grafico:
            st.subheader("Cromatograma (Clica numa letra para editá-ar)")
            st.info("💡 **Dica:** Usa a barra de rolagem abaixo do gráfico para navegar. **Clica em qualquer letra** por cima dos picos para a corrigires.")
            
            fig = go.Figure()
            valor_maximo = 0

            # 1. Desenhar as ondas
            if tracos:
                for base, dados in tracos.items():
                    dados_escalados = [d * escala_vertical for d in dados]
                    if dados_escalados:
                        valor_maximo = max(valor_maximo, max(dados_escalados))
                    
                    fig.add_trace(go.Scatter(
                        y=dados_escalados,
                        name=base,
                        mode='lines',
                        line=dict(color=cores[base], width=1),
                        hoverinfo='skip'
                    ))

            # 2. Desenhar as bases (Clicáveis)
            seq_atual = st.session_state['seq_editada']
            limite = min(len(plocs), len(seq_atual))
            
            # Adicionamos 'customdata' para identificar qual o índice (posição) que foi clicado
            indices = list(range(limite))
            
            fig.add_trace(go.Scatter(
                x=list(plocs)[:limite], 
                y=[valor_maximo * 1.05] * limite, 
                text=list(seq_atual)[:limite],
                mode="text",
                textfont=dict(size=14, color="black", weight="bold"),
                name="Bases (Clicáveis)",
                customdata=indices, 
                hovertext=["Clique para editar" for _ in range(limite)],
                hoverinfo="text"
            ))

            fig.update_layout(
                height=450,
                clickmode='event+select', # Permite seleção por clique
                showlegend=True,
                plot_bgcolor='white',
                margin=dict(l=10, r=10, t=30, b=10),
                xaxis=dict(
                    title="Posição do Traço",
                    rangeslider=dict(visible=True),
                    showgrid=True,
                    gridcolor='lightgrey'
                ),
                yaxis=dict(
                    title="Intensidade",
                    showgrid=True,
                    gridcolor='lightgrey',
                    fixedrange=False
                )
            )

            # 3. Renderizar o gráfico com eventos de seleção ativados (REQUER STREAMLIT >= 1.35)
            selecao = st.plotly_chart(
                fig, 
                use_container_width=True, 
                on_select="rerun",           # O script reinicia quando algo é selecionado
                selection_mode="points",     # Selecionar pontos únicos (as letras)
                key="grafico_sanger"
            )

            # 4. Lógica de Interceção do Clique
            pontos_clicados = selecao.selection.get("points", [])
            # Procurar se algum ponto clicado tem 'customdata' (que são as nossas letras)
            ponto_base = next((p for p in pontos_clicados if "customdata" in p), None)

            if ponto_base:
                # Extrair o índice exato da sequência que o utilizador clicou
                idx_clicado = ponto_base["customdata"]
                base_antiga = seq_atual[idx_clicado]
                
                # Interface de edição cirúrgica!
                st.warning(f"👉 **Modo de Edição:** Selecionaste a base na posição **{idx_clicado + 1}** (Atual: **{base_antiga}**)")
                
                col_input, col_espaco = st.columns([1, 3])
                with col_input:
                    # Campo de texto minúsculo focado apenas numa letra
                    nova_base = st.text_input(
                        "Escreve a nova letra e prime Enter:", 
                        value=base_antiga, 
                        max_chars=1, 
                        key=f"input_base_{idx_clicado}"
                    ).upper()
                    
                    # Se o utilizador alterar a letra, aplicamos na sequência e atualizamos o estado
                    if nova_base and nova_base != base_antiga and nova_base in ['A', 'C', 'T', 'G', 'N']:
                        seq_lista = list(st.session_state['seq_editada'])
                        seq_lista[idx_clicado] = nova_base
                        st.session_state['seq_editada'] = "".join(seq_lista)
                        st.rerun() # Atualiza o gráfico instantaneamente

            # Bloco opcional (para visualização do todo, agora fechado num expansor)
            with st.expander("Ver sequência completa em modo texto"):
                nova_seq_texto = st.text_area(
                    "Sequência Extraída",
                    value=st.session_state['seq_editada'],
                    height=100
                ).upper().strip()
                if nova_seq_texto != st.session_state['seq_editada']:
                    st.session_state['seq_editada'] = nova_seq_texto
                    st.rerun()

        # ==========================================
        # SEPARADOR 2: ALINHAMENTO
        # ==========================================
        with tab_alinhamento:
            st.subheader("Alinhamento com Sequência Teórica")
            
            if seq_referencia:
                if st.button("Executar Alinhamento", type="primary"):
                    aligner = PairwiseAligner()
                    aligner.mode = 'local'
                    aligner.match_score = 2
                    aligner.mismatch_score = -1
                    aligner.open_gap_score = -2
                    aligner.extend_gap_score = -1
                    
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
    st.info("👈 Começa por carregar um ficheiro .ab1 na barra lateral.")
