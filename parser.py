"""Interpretación del enunciado.

1) `analizar_enunciado`: extracción por reglas (sin red) de componentes, cantidades con unidades,
   porcentajes, temperaturas, presiones, reacciones, incógnitas y palabras clave del tipo de proceso.
2) `spec_desde_enunciado`: (opcional) usa un modelo de lenguaje para convertir el enunciado completo
   en la especificación JSON que entiende el solucionador.  Requiere la variable ANTHROPIC_API_KEY.
   El resultado SIEMPRE se valida con el solucionador (grados de libertad, balances) antes de mostrarse.
"""
import json
import os
import re
import urllib.request

from .data import NAMES, molar_mass

_NUM = r"(\d+(?:[.,]\d+)?)"
RE_FLOW = re.compile(_NUM + r"\s*(kg|g|lb|ton|kmol|lbmol|mol)\s*/\s*(s|min|h|hr|d)\b", re.I)
RE_QTY = re.compile(_NUM + r"\s*(kg|g|lb|kmol|lbmol|mol|L|m3|m³)\b(?!\s*/)", re.I)
RE_TEMP = re.compile(_NUM + r"\s*(°\s*[CF]|K\b|grados\s+(?:Celsius|centígrados))", re.I)
RE_PRES = re.compile(_NUM + r"\s*(atm|kPa|bar|mmHg|psi|torr)\b", re.I)
RE_PCT1 = re.compile(_NUM + r"\s*(?:%|por\s+ciento)\s*(?:\(?(?:en\s+)?(masa|peso|mol(?:ar|es)?|volumen)\)?)?\s*(?:de\s+|en\s+)?([A-Za-zÁ-ú][\wÁ-ú]*(?:\s+de\s+carbono)?)", re.I)
RE_PCT2 = re.compile(r"([A-Za-zÁ-ú][\wÁ-ú]*)\s+(?:al|a)\s+" + _NUM + r"\s*%", re.I)
RE_RXN = re.compile(r"((?:\d*\.?\d*\s*[A-Z][A-Za-z0-9()]*\s*\+?\s*)+)\s*(?:->|→|=>)\s*((?:\d*\.?\d*\s*[A-Z][A-Za-z0-9()]*\s*\+?\s*)+)")
RE_FORM = re.compile(r"\b((?:[A-Z][a-z]?\d*){1,8})\b")
QUESTION = re.compile(r"(calcul[ae]\w*|determin[ae]\w*|encuentr[ae]\w*|halle\w*|cu[aá]nt[oa]s?\w*|cu[aá]l(?:es)?|"
                      r"estime\w*|obtenga\w*|what|calculate|find|determine|how much)", re.I)
KEYWORDS = {
    "reacción química": r"reacci[oó]n|reacciona|conversi[oó]n|rendimiento|combusti[oó]n|se quema|exceso|reactor",
    "recirculación (reciclo)": r"recicl|recircul|recycle",
    "purga": r"purga|purge",
    "derivación (bypass)": r"deriv|bypass|by-pass|desv[ií]o",
    "destilación/separación": r"destila|separad|columna|absorb|extracci|cristaliz|secad|evapora",
    "mezclado": r"mezcl|mixer|se combinan",
    "balance de energía": r"calor|entalp|energ[ií]a|temperatura|°\s*[CF]|vapor|enfr[ií]|calent|trabajo",
    "cambio de fase": r"vaporiz|evapora|condens|fusi[oó]n|ebullici|latente|saturad",
    "estado no estacionario/batch": r"por lotes|batch|acumula|semicontinu",
}


def _f(x): return float(x.replace(",", "."))


def analizar_enunciado(texto):
    """Devuelve un dict con lo identificado y un resumen en Markdown (clave 'markdown')."""
    low = texto.lower()
    comps = {}
    for name, formula in NAMES.items():
        if re.search(r"\b" + re.escape(name) + r"\b", low):
            comps[name] = formula
    for m in RE_FORM.finditer(texto):
        tok = m.group(1)
        if any(ch.isdigit() for ch in tok) or tok in ("CO", "NO"):
            try:
                comps.setdefault(tok, tok)
                molar_mass(tok)
            except ValueError:
                comps.pop(tok, None)
    mw = {}
    for n, f in comps.items():
        try:
            mw[n] = molar_mass(f)[0]
        except ValueError:
            pass
    flows = [(_f(a), u.lower(), t.lower()) for a, u, t in RE_FLOW.findall(texto)]
    qtys = [(_f(a), u) for a, u in RE_QTY.findall(texto)]
    temps = [(_f(a), re.sub(r"\s", "", u)) for a, u in RE_TEMP.findall(texto)]
    pres = [(_f(a), u) for a, u in RE_PRES.findall(texto)]
    pcts = []
    for a, base, comp in RE_PCT1.findall(texto):
        pcts.append((comp, _f(a), (base or "no especificada").lower()))
    for comp, a in RE_PCT2.findall(texto):
        pcts.append((comp, _f(a), "no especificada"))
    rxns = [f"{l.strip()} → {r.strip()}" for l, r in RE_RXN.findall(texto)]
    sentences = re.split(r"(?<=[.?!])\s+", texto.strip())
    unknowns = [s.strip() for s in sentences if QUESTION.search(s)]
    conds = [k for k, rx in KEYWORDS.items() if re.search(rx, low)]
    res = dict(componentes=comps, masas_molares=mw, flujos=flows, cantidades=qtys, temperaturas=temps,
               presiones=pres, porcentajes=pcts, reacciones=rxns, incognitas=unknowns, condiciones=conds)
    res["markdown"] = _md(res)
    return res


def _md(r):
    L = ["### Interpretación automática del enunciado (por reglas)\n"]
    L.append("**Componentes:** " + (", ".join(f"{n} ({f}, M={r['masas_molares'].get(n, float('nan')):.3f} g/mol)"
                                              for n, f in r["componentes"].items()) or "no identificados"))
    L.append("**Flujos con unidades:** " + ("; ".join(f"{v:g} {u}/{t}" for v, u, t in r["flujos"]) or "ninguno"))
    L.append("**Cantidades:** " + ("; ".join(f"{v:g} {u}" for v, u in r["cantidades"]) or "ninguna"))
    L.append("**Composiciones (%):** " + ("; ".join(f"{c}: {v:g} % ({b})" for c, v, b in r["porcentajes"]) or "ninguna"))
    L.append("**Temperaturas:** " + ("; ".join(f"{v:g} {u}" for v, u in r["temperaturas"]) or "ninguna"))
    L.append("**Presiones:** " + ("; ".join(f"{v:g} {u}" for v, u in r["presiones"]) or "ninguna"))
    L.append("**Reacciones detectadas:** " + ("; ".join(r["reacciones"]) or "ninguna"))
    L.append("**Condiciones/tipo de proceso:** " + (", ".join(r["condiciones"]) or "no identificadas"))
    L.append("**Incógnitas (lo que se pide):** " + (" | ".join(r["incognitas"]) or "no identificadas"))
    L.append("\n*Esta extracción es heurística: confirme los datos y complete la especificación (corrientes y unidades) "
             "en la pestaña JSON, o use la interpretación con IA.*")
    return "\n\n".join(L)


SCHEMA_PROMPT = """Eres un asistente de Ingeniería Química. Convierte el enunciado de un problema de balance de materia (y opcionalmente energía)
en UN SOLO objeto JSON (sin texto adicional, sin ```), con este esquema:
{
 "title": str, "statement": str (el enunciado), "flow_unit": "kg/h"|"mol/s"|... ,
 "components": [{"name": str, "formula": str, "mw": number(opcional)}],
 "streams": {"<id>": {"name": str, "components": [nombres] (opcional, si se omite son todos),
      "total": "100 kg/h" (opcional), "flows": {"comp": "10 mol/h"} (opcional),
      "fractions": {"comp": 0.4} (opcional), "fraction_basis": "mass"|"mole" (opcional)}},
 "units": [{"name": str, "type": "mixer|splitter|separator|reactor|generic", "inputs": [ids], "outputs": [ids]}],
 "reactions": [{"id": "R1", "unit": nombre, "equation": "CH4 + 2 O2 -> CO2 + 2 H2O",
                "conversion": {"reactant": comp, "value": 0-1, "of_streams": [ids] (opcional)} | "extent": "10 mol/h"}],
 "ratios": [{"a": {"stream": id, "comp": comp|null}, "b": {...}, "value": number, "text": str}],
 "splits": [{"unit": nombre, "stream": id, "fraction": 0-1}],
 "linear": [{"terms": [{"stream": id, "comp": comp|null, "coef": number}], "rhs": number, "text": str}],
 "basis": {"stream": id, "total": "100 mol/h"} (si se elige base de cálculo),
 "energy": {"Tref": 25, "streams": {"id": {"T": °C, "phases": {"comp": "liquid"|"gas"}}}, "units": {"nombre": {"Ws": 0}},
            "thermo": {"formula": {"cp": {"liquid": [a,b,c,d]}, "ref_phase": "liquid", "Hf": number}} (solo si faltan datos)}
}
Reglas: las incógnitas NO llevan datos; cada corriente es entrada de a lo sumo una unidad y salida de a lo sumo una unidad;
todas las relaciones deben ser lineales; use reaction.conversion para conversiones y ratios para excesos de reactivo o razones.
Responde SOLO con el JSON."""


def spec_desde_enunciado(texto, api_key=None, model=None, timeout=90):
    """Convierte el enunciado en especificación JSON con la API de Anthropic (opcional)."""
    api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError("Falta ANTHROPIC_API_KEY (variable de entorno o campo de la interfaz).")
    body = json.dumps({"model": model or os.environ.get("BALANCE_LLM_MODEL", "claude-sonnet-5-5"), "max_tokens": 4000,
                       "system": SCHEMA_PROMPT, "messages": [{"role": "user", "content": texto}]}).encode()
    req = urllib.request.Request("https://api.anthropic.com/v1/messages", data=body, method="POST",
                                 headers={"content-type": "application/json", "x-api-key": api_key,
                                          "anthropic-version": "2023-06-01"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = json.loads(r.read())
    except Exception as e:  # red, clave inválida, etc.
        raise RuntimeError(f"No se pudo consultar la API: {e}")
    txt = "".join(b.get("text", "") for b in data.get("content", []))
    txt = re.sub(r"^```(?:json)?|```$", "", txt.strip(), flags=re.M).strip()
    try:
        return json.loads(txt)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"El modelo no devolvió JSON válido: {e}")
