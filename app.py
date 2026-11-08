"""Interfaz gráfica (Streamlit):  streamlit run app.py"""
import glob
import json
import os

import streamlit as st

from balance_solver import solve_problem
from balance_solver.diagram import draw
from balance_solver.parser import analizar_enunciado, spec_desde_enunciado

st.set_page_config(page_title="Balances de Masa y Energía", page_icon="⚗️", layout="wide")
st.title("⚗️ Balances de Masa y Energía — herramienta educativa")
st.caption("Metodología de Felder & Rousseau: diagrama, grados de libertad, balances, resolución paso a paso y verificación.")

EJ = {os.path.basename(f): f for f in sorted(glob.glob(os.path.join(os.path.dirname(__file__), "examples", "*.json")))}

if "spec_text" not in st.session_state:
    st.session_state.spec_text = ""
if "enunciado" not in st.session_state:
    st.session_state.enunciado = ""

with st.sidebar:
    st.header("Ejemplos resueltos")
    sel = st.selectbox("Cargar un ejemplo", ["—"] + list(EJ))
    if st.button("Cargar ejemplo") and sel != "—":
        with open(EJ[sel], encoding="utf-8") as f:
            data = json.load(f)
        st.session_state.spec_text = json.dumps(data, indent=2, ensure_ascii=False)
        st.session_state.enunciado = data.get("statement", "")
    st.divider()
    st.header("Interpretación con IA (opcional)")
    api_key = st.text_input("ANTHROPIC_API_KEY", type="password", help="Solo se usa para convertir el enunciado en JSON.")
    st.caption("Sin clave, use la extracción por reglas y complete/edite el JSON.")

tab1, tab2, tab3 = st.tabs(["1 · Enunciado", "2 · Especificación (JSON)", "3 · Solución"])

with tab1:
    st.session_state.enunciado = st.text_area("Escriba o pegue el enunciado del ejercicio", st.session_state.enunciado, height=200)
    c1, c2 = st.columns(2)
    if c1.button("🔎 Identificar datos (reglas)"):
        st.markdown(analizar_enunciado(st.session_state.enunciado)["markdown"])
    if c2.button("🤖 Convertir a especificación con IA"):
        try:
            with st.spinner("Interpretando…"):
                spec = spec_desde_enunciado(st.session_state.enunciado, api_key or None)
            st.session_state.spec_text = json.dumps(spec, indent=2, ensure_ascii=False)
            st.success("Especificación generada. Revísela en la pestaña 2 antes de resolver.")
        except Exception as e:
            st.error(str(e))

with tab2:
    st.session_state.spec_text = st.text_area("Especificación del problema (JSON)", st.session_state.spec_text, height=450)
    st.caption("Consulte el README para el formato; los ejemplos de la barra lateral sirven de plantilla.")

with tab3:
    if st.button("▶️ Resolver", type="primary"):
        try:
            spec = json.loads(st.session_state.spec_text)
        except json.JSONDecodeError as e:
            st.error(f"JSON inválido: {e}")
            st.stop()
        try:
            res = solve_problem(spec)
        except (ValueError, KeyError) as e:
            st.error(f"No se pudo plantear el problema: {e}")
            st.stop()
        if res.issues:
            for i in res.issues:
                st.warning(i)
        else:
            st.success("Problema resuelto y verificado.")
        fig = draw(res.problem, res.values if all(s in res.values for s in res.problem.reg.syms) else None, res.energy)
        st.pyplot(fig)
        st.markdown(res.markdown)
        st.download_button("⬇️ Descargar informe (.md)", res.markdown, file_name="informe.md")
