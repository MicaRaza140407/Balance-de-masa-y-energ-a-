"""Utilidades de formato: números y ecuaciones lineales en LaTeX."""
import sympy as sp


def num(x, sig=6):
    x = float(x)
    if x == 0:
        return "0"
    s = f"{x:.{sig}g}"
    if "e" in s:
        m, e = s.split("e")
        return f"{m}\\times10^{{{int(e)}}}"
    return s


def lin_tex(expr, reg, known=None, unit="", rhs="0"):
    """Escribe una expresión lineal (expr = rhs) en LaTeX.

    Las variables presentes en `known` se reemplazan por (valor unidad)."""
    known = known or {}
    expr = sp.expand(expr)
    parts = []
    terms = sorted(sp.Add.make_args(expr), key=lambda t: reg.order_of(t))
    for t in terms:
        coef, rest = t.as_coeff_Mul()
        c = float(coef)
        if rest == 1:
            body = f"{num(abs(c))}\\ \\text{{{unit}}}" if unit else num(abs(c))
        else:
            v = rest
            if v in known:
                val = f"{num(known[v])}\\ \\text{{{unit}}}" if unit else num(known[v])
                body = f"({val})" if abs(c) == 1 else f"{num(abs(c))}({val})"
            else:
                body = reg.tex(v) if abs(c) == 1 else f"{num(abs(c))}\\,{reg.tex(v)}"
        sign = "-" if c < 0 else "+"
        parts.append((sign, body))
    if not parts:
        return f"0 = {rhs}"
    s = ("-" if parts[0][0] == "-" else "") + parts[0][1]
    for sign, body in parts[1:]:
        s += f" {sign} {body}"
    return f"{s} = {rhs}"
