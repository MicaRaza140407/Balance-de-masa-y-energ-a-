"""Masas molares (a partir de la fórmula) y propiedades termodinámicas.

IMPORTANTE: los valores de la tabla THERMO provienen de las tablas B.1 y B.2 del apéndice de
Felder & Rousseau.  Verifíquelos contra su edición del libro; cualquier especie puede
sobrescribirse o añadirse en la sección "thermo" del JSON del problema.
Cp [kJ/(mol·°C)] = a + b·T + c·T² + d·T³  con T en °C.
"""
import re

ATOMIC = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "Cl": 35.45,
          "Na": 22.99, "Ca": 40.08, "P": 30.974, "F": 18.998, "K": 39.098, "Fe": 55.845,
          "Mg": 24.305, "Al": 26.982, "Si": 28.086}

NAMES = {"agua": "H2O", "water": "H2O", "vapor": "H2O", "etanol": "C2H5OH", "ethanol": "C2H5OH",
         "metanol": "CH3OH", "methanol": "CH3OH", "metano": "CH4", "methane": "CH4",
         "propano": "C3H8", "propane": "C3H8", "butano": "C4H10", "etano": "C2H6",
         "oxígeno": "O2", "oxigeno": "O2", "oxygen": "O2", "nitrógeno": "N2", "nitrogeno": "N2",
         "nitrogen": "N2", "hidrógeno": "H2", "hidrogeno": "H2", "hydrogen": "H2",
         "dióxido de carbono": "CO2", "dioxido de carbono": "CO2", "carbon dioxide": "CO2",
         "monóxido de carbono": "CO", "monoxido de carbono": "CO", "amoníaco": "NH3",
         "amoniaco": "NH3", "ammonia": "NH3", "acetona": "C3H6O", "benceno": "C6H6",
         "benzene": "C6H6", "tolueno": "C7H8", "toluene": "C7H8", "glucosa": "C6H12O6"}

THERMO = {
    "H2O": {"cp": {"liquid": [75.4e-3, 0, 0, 0], "gas": [33.46e-3, 0.6880e-5, 0.7604e-8, -3.593e-12]},
            "Tb": 100.0, "Hv": 40.656, "Tm": 0.0, "Hm": 6.01, "ref_phase": "liquid", "Hf": -285.84},
    "N2": {"cp": {"gas": [29.00e-3, 0.2199e-5, 0.5723e-8, -2.871e-12]}, "ref_phase": "gas", "Hf": 0.0},
    "O2": {"cp": {"gas": [29.10e-3, 1.158e-5, -0.6076e-8, 1.311e-12]}, "ref_phase": "gas", "Hf": 0.0},
    "CO2": {"cp": {"gas": [36.11e-3, 4.233e-5, -2.887e-8, 7.464e-12]}, "ref_phase": "gas", "Hf": -393.5},
    "CH4": {"cp": {"gas": [34.31e-3, 5.469e-5, 0.3661e-8, -11.00e-12]}, "ref_phase": "gas", "Hf": -74.85},
}


def parse_formula(formula):
    """Devuelve {elemento: átomos}.  Admite paréntesis: Ca(OH)2."""
    tokens = re.findall(r"([A-Z][a-z]?|\(|\)|\d+)", formula)
    if "".join(tokens) != formula:
        raise ValueError(f"Fórmula no válida: '{formula}'")
    stack = [{}]
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t == "(":
            stack.append({})
        elif t == ")":
            if len(stack) < 2:
                raise ValueError(f"Fórmula no válida: '{formula}'")
            grp = stack.pop()
            mult = 1
            if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                mult = int(tokens[i + 1]); i += 1
            for el, n in grp.items():
                stack[-1][el] = stack[-1].get(el, 0) + n * mult
        elif t.isdigit():
            raise ValueError(f"Fórmula no válida: '{formula}'")
        else:
            n = 1
            if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                n = int(tokens[i + 1]); i += 1
            if t not in ATOMIC:
                raise ValueError(f"Elemento sin masa atómica registrada: {t}")
            stack[-1][t] = stack[-1].get(t, 0) + n
        i += 1
    if len(stack) != 1:
        raise ValueError(f"Fórmula no válida: '{formula}'")
    return stack[0]


def molar_mass(formula):
    """(M en g/mol, texto con el cálculo)."""
    el = parse_formula(formula)
    mw = sum(ATOMIC[e] * n for e, n in el.items())
    expr = " + ".join(f"{n}({ATOMIC[e]})" for e, n in el.items())
    return mw, f"{expr} = {mw:.3f} g/mol"
