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
        
        # Inicialização do estado da sessão
        if 'seq_editada' not in st.session_state or st.session_state.get('id_ficheiro') != ficheiro_carregado.name:
            st.session_state['seq_editada'] = str(record.seq)
            st.session_state['id_ficheiro'] = ficheiro_carregado.name
            st.session_state['cursor_pos'] = 1  # Iniciar o cursor na posição 1

        tab_grafico, tab_alinhamento = st.tabs(["📊 Cromatograma Interativo", "🔍 Alinhamento Global"])

        with tab_grafico:
            st.subheader("Modo de Edição (Estilo SnapGene)")
            st.info("💡 Usa os controlos abaixo para mover o **Cursor Vermelho** no gráfico. Assim que estiver no pico desejado, altera a base na caixa.")
            
            seq_atual = st.session_state['seq_editada']
            limite = min(len(plocs), len(seq_atual)) if plocs else len(seq_atual)
            
            # --- PAINEL DE CONTROLO DO CURSOR ---
            st.markdown("---")
            col_nav1, col_nav2, col_edit, col_vazio = st.columns([1.5, 2, 2, 2])
            
            with col_nav1:
                st.write("**Navegação Rápida**")
                btn_esq, btn_dir = st.columns(2)
                if btn_esq.button("⬅️ Ant."):
                    st.session_state['cursor_pos'] = max(1, st.session_state['cursor_pos'] - 1)
                if btn_dir.button("Seg. ➡️"):
                    st.session_state['cursor_pos'] = min(limite, st.session_state['cursor_pos'] + 1)
            
            with col_nav2:
                # O widget number_input atualiza diretamente a variável de sessão 'cursor_pos'
                st.number_input(
                    "📍 Ir para a Posição:", 
                    min_value=1, 
                    max_value=limite, 
                    key="cursor_pos"
                )
            
            with col_edit:
                # Posição atual baseada no cursor (1-indexado para o utilizador, 0-indexado para a lista)
                idx_atual = st.session_state['cursor_pos'] - 1
                base_atual = seq_atual[idx_atual]
                
                nova_base = st.text_input(
                    f"Substituir base {idx_atual + 1}:", 
                    value=base_atual, 
                    max_chars=1
                ).upper()
                
                # Guardar a edição
                if nova_base and nova_base != base_atual and nova_base in ['A', 'C', 'T', 'G', 'N']:
                    seq_lista = list(st.session_state['seq_editada'])
                    seq_lista[idx_atual] = nova_base
                    st.session_state['seq_editada'] = "".join(seq_lista)
                    st.rerun()
            st.markdown("---")

            # --- CONSTRUÇÃO DO GRÁFICO PLOTLY ---
            fig = go.Figure()
            valor_maximo = 0

            # 1. Desenhar as ondas (traces)
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

            # 2. Desenhar todas as bases no topo dos picos
            fig.add_trace(go.Scatter(
                x=list(plocs)[:limite], 
                y=[valor_maximo * 1.05] * limite, 
                text=list(seq_atual)[:limite],
                mode="text",
                textfont=dict(size=14, color="black"),
                name="Bases Chamadas"
            ))

            # 3. Desenhar o CURSOR VERMELHO (Estilo SnapGene)
            if plocs and idx_atual < len(plocs):
                x_cursor = plocs[idx_atual]
                
                # Linha vertical
                fig.add_vline(x=x_cursor, line_width=2, line_dash="dash", line_color="rgba(255, 0, 0, 0.7)")
                
                # Destaque na letra atual (Caixa vermelha)
                fig.add_annotation(
                    x=x_cursor,
                    y=valor_maximo * 1.05,
                    text=seq_atual[idx_atual],
                    showarrow=False,
                    font=dict(color="white", size=16, weight="bold"),
                    bgcolor="red",
                    bordercolor="darkred",
                    borderwidth=2,
                    borderpad=4
                )

            # 4. Configurações de layout (Zoom e Pan)
            fig.update_layout(
                height=450,
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

            st.plotly_chart(fig, use_container_width=True)

            # --- EDITOR EM MASSA ---
            with st.expander("🛠️ Ver/Editar a sequência completa em modo de texto"):
                nova_seq_texto = st.text_area(
                    "Podes colar uma sequência inteira aqui se preferires:",
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
