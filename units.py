"""Unidades y conversiones (masa, moles, tiempo).

Todas las conversiones pasan por una unidad base: kg para masa, mol para moles.
Convertir entre masa y moles requiere la masa molar (g/mol).
"""
import re

MASS_U = {"kg": 1.0, "g": 1e-3, "lb": 0.45359237, "ton": 1000.0}
MOLE_U = {"mol": 1.0, "kmol": 1000.0, "lbmol": 453.59237}
SEC = {"s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0}
# Unidad molar equivalente cuando hay que pasar de base másica a molar
MOLE_OF = {"kg": "kmol", "g": "mol", "lb": "lbmol", "ton": "kmol",
           "mol": "mol", "kmol": "kmol", "lbmol": "lbmol"}
_ALIASES = {"hr": "h", "hora": "h", "horas": "h", "seg": "s", "sec": "s", "dia": "d", "día": "d",
            "t": "ton", "lbs": "lb", "gmol": "mol", "g-mol": "mol", "kgmol": "kmol", "lb-mol": "lbmol"}


def norm(u):
    u = u.strip()
    return _ALIASES.get(u.lower(), u.lower())


def kind_of(amount_unit):
    if amount_unit in MASS_U:
        return "mass"
    if amount_unit in MOLE_U:
        return "mole"
    raise ValueError(f"Unidad de cantidad desconocida: '{amount_unit}'. Use kg, g, lb, ton, mol, kmol o lbmol.")


def base_factor(u):
    return MASS_U[u] if u in MASS_U else MOLE_U[u]


def amount_factor(from_u, to_u, mw=None):
    """Factor f tal que (cantidad en from_u) * f = (cantidad en to_u).  mw en g/mol."""
    kf, kt = kind_of(from_u), kind_of(to_u)
    if kf == kt:
        return base_factor(from_u) / base_factor(to_u)
    if mw is None:
        raise ValueError("Se necesita la masa molar para convertir entre masa y moles.")
    if kf == "mass":      # kg -> mol : x1000/MW
        return base_factor(from_u) * 1000.0 / mw / base_factor(to_u)
    return base_factor(from_u) * mw / 1000.0 / base_factor(to_u)   # mol -> kg


def time_factor(from_t, to_t):
    """Factor para convertir un flujo por from_t en un flujo por to_t."""
    return SEC[to_t] / SEC[from_t]


_QTY = re.compile(r"^\s*([-+]?\d*[.,]?\d+(?:[eE][-+]?\d+)?)\s*([A-Za-z\-]+)\s*(?:/\s*([A-Za-z]+))?\s*$")


def parse_quantity(q):
    """Devuelve (valor, unidad_cantidad|None, unidad_tiempo|None).

    Acepta número, "100 kg/h", "5 kmol" o {"value": 100, "unit": "kg/h"}."""
    if isinstance(q, (int, float)):
        return float(q), None, None
    if isinstance(q, dict):
        if q.get("unit"):
            return parse_quantity(f"{q['value']} {q['unit']}")
        return float(q["value"]), None, None
    if isinstance(q, str):
        s = q.strip()
        try:
            return float(s.replace(",", ".")), None, None
        except ValueError:
            pass
        m = _QTY.match(s)
        if m:
            a = norm(m.group(2))
            t = norm(m.group(3)) if m.group(3) else None
            kind_of(a)
            if t is not None and t not in SEC:
                raise ValueError(f"Unidad de tiempo desconocida en '{q}'.")
            return float(m.group(1).replace(",", ".")), a, t
    raise ValueError(f"No se pudo interpretar la cantidad: {q!r}")


def parse_flow_unit(u):
    a, _, t = u.partition("/")
    a, t = norm(a), norm(t)
    kind_of(a)
    if t not in SEC:
        raise ValueError(f"La unidad de flujo '{u}' debe tener la forma cantidad/tiempo (p. ej. kg/h).")
    return a, t
