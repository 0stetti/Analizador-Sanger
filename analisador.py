import streamlit as st
from Bio import SeqIO
from Bio.Align import PairwiseAligner, AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import plotly.graph_objects as go
import numpy as np
import io

# Configuração da página
st.set_page_config(page_title="Sanger Editor & Aligner", layout="wide")

# --- GERENCIAMENTO DE ESTADO (SESSION STATE) ---
# Isso é crucial para manter as edições do usuário na memória
if 'reads_data' not in st.session_state:
    st.session_state['reads_data'] = {} # Dicionário para guardar {nome_arquivo: sequencia_editada}

def get_chromatogram_traces(record):
    """Extrai os canais de cor do arquivo AB1"""
    channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12'] # G, A, T, C (Padrão ABI)
    colors = {'G': 'black', 'A': 'green', 'T': 'red', 'C': 'blue'}
    bases = ['G', 'A', 'T', 'C']
    traces = {}
    
    # Tenta extrair os canais raw do arquivo
    if 'abif_raw' in record.annotations:
        raw = record.annotations['abif_raw']
        for i, channel in enumerate(channels):
            if channel in raw:
                traces[bases[i]] = list(raw[channel])
    return traces

def plot_chromatogram(traces, sequence_length):
    """Gera o gráfico interativo com Plotly"""
    fig = go.Figure()
    
    # Adiciona cada canal ao gráfico
    colors = {'G': 'black', 'A': 'green', 'T': 'red', 'C': 'blue'}
    for base, data in traces.items():
        if data:
            fig.add_trace(go.Scatter(
                y=data, 
                mode='lines', 
                name=base, 
                line=dict(color=colors[base], width=1),
                opacity=0.7
            ))
            
    fig.update_layout(
        title="Visualização do Cromatograma (Zoom Interativo)",
        height=300,
        margin=dict(l=20, r=20, t=40, b=20),
        xaxis_title="Posição do Trace",
        yaxis_title="Intensidade de Fluorescência",
        legend_title="Bases"
    )
    return fig

# --- INTERFACE PRINCIPAL ---

st.title("🧬 Sanger Master: Edição e Alinhamento")

# 1. INPUTS NA SIDEBAR
with st.sidebar:
    st.header("1. Upload e Referência")
    
    # Upload Múltiplo
    uploaded_files = st.file_uploader("Subir arquivos .ab1", type=["ab1"], accept_multiple_files=True)
    
    # Processar novos arquivos
    if uploaded_files:
        for uploaded_file in uploaded_files:
            if uploaded_file.name not in st.session_state['reads_data']:
                try:
                    # Ler o arquivo binário
                    bytes_data = uploaded_file.read()
                    record = SeqIO.read(io.BytesIO(bytes_data), "abi")
                    
                    # Salvar no estado
                    st.session_state['reads_data'][uploaded_file.name] = {
                        'original_seq': str(record.seq),
                        'edited_seq': str(record.seq), # Começa igual a original
                        'traces': get_chromatogram_traces(record),
                        'record_obj': record
                    }
                except Exception as e:
                    st.error(f"Erro ao ler {uploaded_file.name}: {e}")
            # Resetar ponteiro do arquivo para não dar erro se o streamlit recarregar
            uploaded_file.seek(0)

    st.markdown("---")
    ref_seq_input = st.text_area("Sequência Teórica (Referência)", height=150, placeholder="Cole aqui a sequência do gene esperado (ex: mclover 3)...").upper().strip()

# --- ABAS DA APLICAÇÃO ---
tab_editor, tab_results = st.tabs(["✏️ Editor de Sequências", "📊 Alinhamento e Consenso"])

with tab_editor:
    st.subheader("Visualizar e Editar Bases")
    
    if not st.session_state['reads_data']:
        st.info("Faça upload de arquivos .ab1 na barra lateral para começar.")
    
    # Loop pelos arquivos carregados para criar editores individuais
    for filename, data in st.session_state['reads_data'].items():
        with st.expander(f"Arquivo: {filename}", expanded=True):
            
            # 1. Plota o gráfico (Cromatograma)
            if data['traces']:
                st.plotly_chart(plot_chromatogram(data['traces'], len(data['original_seq'])), use_container_width=True)
            else:
                st.warning("Dados de traço (traces) não encontrados neste arquivo.")

            # 2. Área de Edição
            col1, col2 = st.columns([3, 1])
            with col1:
                # Onde a mágica acontece: O usuário edita, e nós salvamos no Session State
                new_seq = st.text_area(
                    f"Sequência de Bases ({filename})", 
                    value=data['edited_seq'], 
                    height=100,
                    key=f"editor_{filename}" # Chave única para persistência
                ).upper().strip()
                
                # Atualiza o estado se houve mudança
                if new_seq != data['edited_seq']:
                    st.session_state['reads_data'][filename]['edited_seq'] = new_seq
                    st.toast(f"Sequência de {filename} atualizada!", icon="✅")
            
            with col2:
                st.write("Estatísticas:")
                st.caption(f"Tamanho Original: {len(data['original_seq'])} bp")
                st.caption(f"Tamanho Atual: {len(new_seq)} bp")
                if new_seq != data['original_seq']:
                    st.caption("⚠️ **Editado Manualmente**")

with tab_results:
    st.subheader("Alinhamento Global")
    
    if ref_seq_input and st.session_state['reads_data']:
        if st.button("Rodar Alinhamento"):
            aligner = PairwiseAligner()
            aligner.mode = 'local' # Local é melhor para achar a região de interesse
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -1
            
            st.markdown("### Resultados por Amostra")
            
            valid_sequences = [] # Lista para gerar consenso
            
            # Alinha cada sequência editada contra a referência
            for filename, data in st.session_state['reads_data'].items():
                seq_usuario = data['edited_seq']
                alignment = aligner.align(ref_seq_input, seq_usuario)[0]
                
                st.markdown(f"**Amostra:** `{filename}` | Score: {alignment.score}")
                st.code(str(alignment), language='text')
                
                # Adiciona à lista para tentar consenso (somente se tiver tamanho similar)
                # Nota: Consenso real requer Multiple Sequence Alignment (MSA), 
                # aqui faremos uma aproximação visual empilhando os dados
                valid_sequences.append(SeqRecord(Seq(seq_usuario), id=filename))
            
            st.divider()
            st.markdown("### Resumo / Consenso Simplificado")
            st.info("Para um consenso verdadeiro, seria ideal usar algoritmos de MSA (como Muscle). Abaixo, uma visão das sequências carregadas:")
            
            # Mostra as sequências brutas uma em cima da outra para comparação visual rápida
            df_display = {
                "Nome do Arquivo": [],
                "Sequência (Início)": [],
                "Sequência (Fim)": []
            }
            for v in valid_sequences:
                df_display["Nome do Arquivo"].append(v.id)
                df_display["Sequência (Início)"].append(str(v.seq)[:50] + "...")
                df_display["Sequência (Fim)"].append("..." + str(v.seq)[-50:])
            
            st.dataframe(df_display)

    else:
        st.warning("Por favor, certifique-se de ter carregado arquivos .ab1 E inserido uma sequência de referência.")
