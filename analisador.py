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
    Extrai os traços (ondas), a sequência original, a posição dos picos (PLOC) 
    e a qualidade (Phred).
    """
    canais = ['DATA9', 'DATA10', 'DATA11', 'DATA12'] 
    mapa_bases = {'DATA9': 'G', 'DATA10': 'A', 'DATA11': 'T', 'DATA12': 'C'}
    
    # Cores pasteis suaves para as linhas e fundo
    cores_linhas = {'G': '#9E9E9E', 'A': '#81C784', 'T': '#E57373', 'C': '#64B5F6'}
    cores_fundo = {'G': 'rgba(158,158,158,0.2)', 'A': 'rgba(129,199,132,0.2)', 'T': 'rgba(229,115,115,0.2)', 'C': 'rgba(100,181,246,0.2)'}
    
    tracos = {}
    
    if 'abif_raw' in record.annotations:
        raw = record.annotations['abif_raw']
        
        for canal in canais:
            if canal in raw:
                base = mapa_bases[canal]
                tracos[base] = list(raw[canal])
        
        ploc = raw.get('PLOC_1', raw.get('PLOC_2', []))
        
        # Tentar extrair qualidade Phred
        qualidade = record.letter_annotations.get("phred_quality", [0]*len(ploc))
        
        return tracos, ploc, qualidade, cores_linhas, cores_fundo
    return None, None, None, None, None

def cor_da_letra(base):
    # Mantemos as letras com o mesmo tom pastel das linhas
    cores_pasteis = {
        'A': '#81C784', 
        'T': '#E57373', 
        'C': '#64B5F6', 
        'G': '#9E9E9E'
    }
    return cores_pasteis.get(base, 'gray')

# --- INTERFACE PRINCIPAL ---
st.title("🧬 Sanger Pro: Visualização e Alinhamento")

with st.sidebar:
    st.header("1. Carregar Dados")
    ficheiro_carregado = st.file_uploader("Ficheiro .ab1", type=["ab1"])
    
    st.markdown("---")
    st.header("2. Controlos de Visualização")
    escala_vertical = st.slider("Amplitude dos Picos", 1.0, 10.0, 1.0, 0.1)
    zoom_horizontal = st.slider("Bases por vista (Zoom X)", 20, 200, 80, 10)
    
    st.markdown("---")
    st.header("3. Referência (Opcional)")
    seq_referencia = st.text_area("Sequência Teórica (ex: mclover 3)", height=150).upper().strip()

if ficheiro_carregado:
    try:
        dados_bytes = ficheiro_carregado.read()
        record = SeqIO.read(io.BytesIO(dados_bytes), "abi")
        
        tracos, plocs, qualidade, cores, cores_fundo = obter_dados_sanger(record)
        
        # Inicialização do estado da sessão
        if 'seq_editada' not in st.session_state or st.session_state.get('id_ficheiro') != ficheiro_carregado.name:
            st.session_state['seq_editada'] = str(record.seq)
            st.session_state['id_ficheiro'] = ficheiro_carregado.name
            st.session_state['cursor_pos'] = 1  

        tab_grafico, tab_alinhamento = st.tabs(["📊 Cromatograma", "🔍 Alinhamento"])

        with tab_grafico:
            st.subheader("Modo de Edição (Estilo SnapGene)")
            
            seq_atual = st.session_state['seq_editada']
            limite = min(len(plocs), len(seq_atual)) if plocs else len(seq_atual)
            
            # --- PAINEL DE CONTROLO DO CURSOR ---
            col_nav1, col_nav2, col_edit, col_vazio = st.columns([1.5, 2, 2, 2])
            
            with col_nav1:
                st.write("**Navegação Rápida**")
                btn_esq, btn_dir = st.columns(2)
                if btn_esq.button("⬅️ Ant."):
                    st.session_state['cursor_pos'] = max(1, st.session_state['cursor_pos'] - 1)
                if btn_dir.button("Seg. ➡️"):
                    st.session_state['cursor_pos'] = min(limite, st.session_state['cursor_pos'] + 1)
            
            with col_nav2:
                st.number_input(
                    "📍 Ir para a Posição:", 
                    min_value=1, 
                    max_value=limite, 
                    key="cursor_pos"
                )
            
            with col_edit:
                idx_atual = st.session_state['cursor_pos'] - 1
                base_atual = seq_atual[idx_atual]
                
                nova_base = st.text_input(
                    f"Substituir base {idx_atual + 1}:", 
                    value=base_atual, 
                    max_chars=1
                ).upper()
                
                if nova_base and nova_base != base_atual and nova_base in ['A', 'C', 'T', 'G', 'N', '-']:
                    seq_lista = list(st.session_state['seq_editada'])
                    seq_lista[idx_atual] = nova_base
                    st.session_state['seq_editada'] = "".join(seq_lista)
                    st.rerun()
            st.markdown("---")

            # --- CONSTRUÇÃO DO GRÁFICO PLOTLY (ESTILO SNAPGENE) ---
            fig = go.Figure()
            valor_maximo = 0

            # 1. Desenhar as ondas com preenchimento (estilo SnapGene)
            if tracos:
                for base, dados in tracos.items():
                    dados_escalados = [d * escala_vertical for d in dados]
                    if dados_escalados:
                        valor_maximo = max(valor_maximo, max(dados_escalados))
                    
                    fig.add_trace(go.Scatter(
                        y=dados_escalados,
                        name=f"Canal {base}",
                        mode='lines',
                        line=dict(color=cores[base], width=1.5),
                        fill='tozeroy',           # Preenchimento
                        fillcolor=cores_fundo[base], # Cor transparente
                        hoverinfo='skip'
                    ))

            # Preparar as letras coloridas
            lista_plocs = list(plocs)[:limite]
            lista_letras = list(seq_atual)[:limite]
            cores_letras = [cor_da_letra(b) for b in lista_letras]
            
            # Hover text com a qualidade (Phred)
            textos_hover = [f"Posição: {i+1}<br>Qualidade (Phred): {qualidade[i] if i < len(qualidade) else 'N/A'}" for i in range(limite)]

            # 2. Desenhar todas as letras no topo
            fig.add_trace(go.Scatter(
                x=lista_plocs, 
                y=[valor_maximo * 1.05] * limite, 
                text=lista_letras,
                mode="text",
                textfont=dict(size=14, color=cores_letras, family="monospace", weight="bold"),
                name="Sequência",
                hovertext=textos_hover,
                hoverinfo="text"
            ))

            # 3. Adicionar marcadores verticais a cada 10 bases (Linha do tempo)
            for i in range(9, limite, 10):
                fig.add_vline(x=plocs[i], line_width=1, line_dash="dot", line_color="rgba(128,128,128,0.5)")
                fig.add_annotation(
                    x=plocs[i], y=0, 
                    text=str(i+1), showarrow=False, 
                    yshift=-20, font=dict(color="gray", size=10)
                )

            # 4. Desenhar o CURSOR (Sombra / Destaque)
            if plocs and idx_atual < len(plocs):
                x_cursor = plocs[idx_atual]
                
                # Fundo azul claro no pico atual (como no SnapGene)
                largura_pico = 15 # estimativa de largura
                fig.add_vrect(
                    x0=x_cursor - largura_pico, x1=x_cursor + largura_pico,
                    fillcolor="rgba(0, 150, 255, 0.2)",
                    layer="below", line_width=0,
                )
                
                # Destacar a letra com uma caixa no topo
                fig.add_annotation(
                    x=x_cursor,
                    y=valor_maximo * 1.05,
                    text=seq_atual[idx_atual],
                    showarrow=False,
                    font=dict(color="white", size=16, weight="bold"),
                    bgcolor="rgba(0, 100, 255, 0.8)",
                    bordercolor="darkblue", borderwidth=1, borderpad=3
                )

            # 5. Cálculo Dinâmico da Janela de Zoom (Para não espremer tudo)
            centro_x = plocs[idx_atual] if (plocs and idx_atual < len(plocs)) else 0
            # Mostrar metade das bases pedidas para a esquerda e metade para a direita
            raio_zoom = zoom_horizontal * 15 # 15 é a distância média aproximada entre picos
            
            x_min = max(0, centro_x - raio_zoom)
            x_max = centro_x + raio_zoom

            # 6. Layout Final
            fig.update_layout(
                height=450,
                showlegend=False, # Ocultar legenda para ter mais espaço limpo
                plot_bgcolor='white',
                margin=dict(l=10, r=10, t=30, b=40),
                xaxis=dict(
                    title="Posição do Traço",
                    range=[x_min, x_max], # <- A MÁGICA ESTÁ AQUI (Auto-Zoom no Cursor)
                    rangeslider=dict(visible=True, thickness=0.08),
                    showgrid=False,
                    zeroline=False
                ),
                yaxis=dict(
                    title="Intensidade",
                    showgrid=True,
                    gridcolor='rgba(200,200,200,0.2)',
                    fixedrange=False,
                    zeroline=False
                )
            )

            st.plotly_chart(fig, use_container_width=True)

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
