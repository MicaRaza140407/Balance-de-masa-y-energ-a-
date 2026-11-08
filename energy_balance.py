"""Balances de energía (Felder & Rousseau, caps. 7-9) para sistemas abiertos en estado estacionario.

Ecuación base (ΔEc y ΔEp despreciables):      Q − Ws = ΔH
con ΔH = Σ ξ_r·ΔH°_r(25 °C) + Σ n_sal·Ĥ_sal − Σ n_ent·Ĥ_ent   (método del avance de reacción).
Estados de referencia: cada especie en su fase de referencia a Tref (por defecto 25 °C, 1 atm),
que son los mismos estados en los que se tabulan los calores de formación.
Convención de signos: Q > 0 si el calor entra al sistema; Ws > 0 si el sistema realiza trabajo.
"""
from scipy.optimize import brentq
import copy

from . import units as U
from .data import THERMO
from .render import num


def _poly_int(coefs, T1, T2):
    a, b, c, d = coefs
    return (a * (T2 - T1) + b / 2 * (T2 ** 2 - T1 ** 2) + c / 3 * (T2 ** 3 - T1 ** 3) + d / 4 * (T2 ** 4 - T1 ** 4))


def _poly_tex(coefs, T1, T2, label):
    a, b, c, d = coefs
    return (f"\\int_{{{num(T1)}}}^{{{num(T2)}}} C_{{p,{label}}}\\,dT = "
            f"{num(a)}({num(T2)}-{num(T1)})"
            + (f" + \\tfrac{{{num(b)}}}{{2}}({num(T2)}^2-{num(T1)}^2)" if b else "")
            + (f" + \\tfrac{{{num(c)}}}{{3}}({num(T2)}^3-{num(T1)}^3)" if c else "")
            + (f" + \\tfrac{{{num(d)}}}{{4}}({num(T2)}^4-{num(T1)}^4)" if d else "")
            + f" = {num(_poly_int(coefs, T1, T2))}\\ \\text{{kJ/mol}}")


def enthalpy(name, th, phase, T, Tref, steps=None):
    """Ĥ [kJ/mol] de la especie en (fase, T) respecto a (fase de referencia, Tref)."""
    ref = th["ref_phase"]
    cp = th["cp"]
    def need(ph):
        if ph not in cp:
            raise ValueError(f"Falta Cp en fase '{ph}' para {name}. Agréguelo en la sección 'thermo'.")
        return cp[ph]
    log = steps if steps is not None else []
    if phase == ref:
        v = _poly_int(need(phase), Tref, T)
        log.append(f"{name} ({phase}, sin cambio de fase): $\\hat H = {_poly_tex(need(phase), Tref, T, phase[0])}$")
        return v
    if "Tb" not in th or "Hv" not in th:
        raise ValueError(f"Para {name} con cambio de fase se necesitan 'Tb' y 'Hv'.")
    Tb, Hv = th["Tb"], th["Hv"]
    if {phase, ref} != {"liquid", "gas"}:
        raise ValueError(f"Cambio de fase {ref}→{phase} no soportado para {name} (solo líquido↔gas).")
    sgn = 1 if ref == "liquid" else -1
    h1 = _poly_int(need(ref), Tref, Tb)
    h2 = sgn * Hv
    h3 = _poly_int(need(phase), Tb, T)
    log.append(f"{name} ({ref}→{phase}) se calcula por trayectoria: {ref} de {num(Tref)} °C a T_b = {num(Tb)} °C, "
               f"cambio de fase a T_b y {phase} de T_b a {num(T)} °C.")
    log.append(f"$\\hat H_1 = {_poly_tex(need(ref), Tref, Tb, ref[0])}$")
    log.append(f"$\\hat H_2 = {'+' if sgn > 0 else '-'}\\Delta\\hat H_v({num(Tb)}\\,°C) = {num(h2)}\\ \\text{{kJ/mol}}$")
    log.append(f"$\\hat H_3 = {_poly_tex(need(phase), Tb, T, phase[0])}$")
    log.append(f"$\\hat H = \\hat H_1 + \\hat H_2 + \\hat H_3 = {num(h1)} + {num(h2)} + {num(h3)} = {num(h1 + h2 + h3)}\\ \\text{{kJ/mol}}$")
    return h1 + h2 + h3


def phase_warnings(name, th, phase, T, stream):
    out = []
    if "Tb" in th:
        if phase == "liquid" and T > th["Tb"] + 1e-9:
            out.append(f"Corriente {stream}: {name} como líquido a {T:g} °C supera su punto de ebullición normal "
                       f"({th['Tb']:g} °C a 1 atm): revise la fase o la presión.")
        if phase == "gas" and T < th["Tb"] - 1e-9 and name == "H2O":
            out.append(f"Corriente {stream}: {name} como gas a {T:g} °C está por debajo de su T de ebullición a 1 atm; "
                       f"solo es válido si la presión parcial es menor que la de saturación.")
    return out


def solve_energy(problem, values, steps, issues):
    """Balance de energía por unidad usando los flujos de materia ya resueltos.

    Devuelve dict {unidad: {Q, Ws, dH}} en kJ por unidad de tiempo del problema."""
    spec = problem.spec["energy"]
    Tref = float(spec.get("Tref", 25))
    thermo = copy.deepcopy(THERMO)
    for k, v in spec.get("thermo", {}).items():
        thermo.setdefault(k, {}).update(v)
    th_of = lambda c: thermo.get(problem.comps[c]["formula"]) or thermo.get(c) or (_ for _ in ()).throw(
        ValueError(f"No hay datos termodinámicos de '{c}'. Agréguelos en energy.thermo."))
    tf = problem.time_u
    results = {}
    steps.append("## Balance de energía")
    steps.append("**Método (Felder, caps. 7–9).** Sistema abierto en estado estacionario con ΔEc y ΔEp despreciables: "
                 "$$\\dot Q - \\dot W_s = \\Delta\\dot H$$ "
                 "Las entalpías específicas se calculan respecto a un estado de referencia: cada especie en su fase "
                 f"natural a {Tref:g} °C y 1 atm. Con ese criterio la entalpía de una corriente que sale es "
                 "$\\sum \\dot n_i\\hat H_i$ y, si hay reacción, se agrega el término $\\xi\\,\\Delta H^\\circ_r$ "
                 "(método del avance de reacción, adecuado porque ξ ya se obtuvo en el balance de materia). "
                 "Se desprecian los efectos de presión sobre la entalpía.")
    # calores de reacción
    dHr = {}
    for r in problem.reactions:
        parts, tot = [], 0.0
        for c, nu in r["stoich"].items():
            hf = th_of(c).get("Hf")
            if hf is None:
                raise ValueError(f"Falta ΔHf de {c} (energy.thermo).")
            tot += float(nu) * hf
            parts.append(f"({num(float(nu))})({num(hf)})")
        dHr[r["id"]] = tot
        steps.append(f"**Calor estándar de la reacción {r['id']}** ({r['equation']}), a partir de los calores de formación: "
                     f"$$\\Delta H^\\circ_{{{r['id']}}} = \\sum_i \\nu_i\\,\\Delta\\hat H^\\circ_{{f,i}} = "
                     f"{' + '.join(parts)} = {num(tot)}\\ \\text{{kJ/mol}}$$")

    def n_mol(s, c):
        v = float(values[problem.var[(s, c)]])
        return v * U.amount_factor(problem.int_u, "mol", problem.comps[c]["mw"] if problem.kind == "mass" else None)

    def stream_H(s, label, steps_):
        info = spec["streams"][s]
        T = float(info["T"])
        tot = 0.0
        steps_.append(f"*Corriente {s}* a {T:g} °C:")
        for c in problem.streams[s]["components"]:
            n = n_mol(s, c)
            if abs(n) < 1e-12:
                continue
            th = th_of(c)
            ph = info.get("phases", {}).get(c, th["ref_phase"])
            issues.extend(phase_warnings(c, th, ph, T, s))
            sub = []
            h = enthalpy(c, th, ph, T, Tref, sub)
            steps_.extend(f"- {x}" for x in sub)
            steps_.append(f"- $\\dot n_{{\\mathrm{{{c}}},{s}}}\\hat H = ({num(n)}\\ \\text{{mol/{tf}}})({num(h)}\\ \\text{{kJ/mol}}) = {num(n * h)}\\ \\text{{kJ/{tf}}}$")
            tot += n * h
        steps_.append(f"- $\\dot H_{{{s}}} = {num(tot)}\\ \\text{{kJ/{tf}}}$")
        return tot

    solveT = spec.get("solve_T")
    for uname, uspec in spec.get("units", {}).items():
        ins, outs = problem.unit_streams[uname]
        Ws = float(uspec.get("Ws", 0))
        steps.append(f"### Unidad «{uname}»")
        rxs = [r for r in problem.reactions if r.get("unit") == uname]
        if solveT and solveT["unit"] == uname:
            Qgiven = float(solveT["Q"])
            s_t = str(solveT["stream"])
            lo, hi = solveT.get("bracket", [-100.0, 3000.0])

            def f(T):
                spec["streams"][s_t]["T"] = T
                tot = 0.0
                for s in outs:
                    tot += _quiet_H(problem, s, spec, th_of, Tref, n_mol)
                for s in ins:
                    tot -= _quiet_H(problem, s, spec, th_of, Tref, n_mol)
                for r in rxs:
                    tot += float(values[problem.xi[r["id"]]]) * dHr[r["id"]]
                return tot - (Qgiven - Ws)
            T_sol = brentq(f, lo, hi)
            spec["streams"][s_t]["T"] = T_sol
            steps.append(f"Se conoce Q = {num(Qgiven)} kJ/{tf} y se busca la temperatura de la corriente {s_t}. "
                         f"La ecuación $Q-W_s=\\Delta H$ es no lineal en T (Cp depende de T); se resuelve numéricamente "
                         f"(método de Brent) en [{lo:g}, {hi:g}] °C: $T_{{{s_t}}} = {num(T_sol)}\\ °C$.")
        H_in, H_out = 0.0, 0.0
        steps.append("**Entalpías de las corrientes de entrada**")
        for s in ins:
            H_in += stream_H(s, "ent", steps)
        steps.append("**Entalpías de las corrientes de salida**")
        for s in outs:
            H_out += stream_H(s, "sal", steps)
        react = 0.0
        for r in rxs:
            xi = float(values[problem.xi[r["id"]]]) * U.amount_factor(
                U.MOLE_OF[problem.amount_u] if problem.kind == "mass" else problem.int_u, "mol")
            react += xi * dHr[r["id"]]
            steps.append(f"Término de reacción: $\\xi_{{{r['id']}}}\\Delta H^\\circ_{{{r['id']}}} = ({num(xi)}\\ \\text{{mol/{tf}}})"
                         f"({num(dHr[r['id']])}\\ \\text{{kJ/mol}}) = {num(xi * dHr[r['id']])}\\ \\text{{kJ/{tf}}}$")
        dH = react + H_out - H_in
        Q = dH + Ws
        rtxt = f"{num(react)} + " if rxs else ""
        xi_term = "\\xi\\Delta H_r^\\circ + " if rxs else ""
        steps.append(f"$$\\Delta\\dot H = {xi_term}\\sum \\dot H_{{sal}} - \\sum \\dot H_{{ent}} = "
                     f"{rtxt}({num(H_out)}) - ({num(H_in)}) = {num(dH)}\\ \\text{{kJ/{tf}}}$$")
        steps.append(f"$$\\dot Q = \\Delta\\dot H + \\dot W_s = ({num(dH)}) + ({num(Ws)}) = {num(Q)}\\ \\text{{kJ/{tf}}}"
                     f" = {num(Q / U.SEC[tf])}\\ \\text{{kW}}$$")
        steps.append(f"**Resultado:** el sistema {'absorbe' if Q > 0 else 'libera'} {num(abs(Q))} kJ/{tf} "
                     f"({'calor que debe suministrarse' if Q > 0 else 'calor que debe retirarse'}).")
        results[uname] = dict(Q=Q, Ws=Ws, dH=dH, H_in=H_in, H_out=H_out, react=react)
    return results


def _quiet_H(problem, s, spec, th_of, Tref, n_mol):
    info = spec["streams"][s]
    T = float(info["T"])
    tot = 0.0
    for c in problem.streams[s]["components"]:
        n = n_mol(s, c)
        if abs(n) < 1e-12:
            continue
        th = th_of(c)
        ph = info.get("phases", {}).get(c, th["ref_phase"])
        tot += n * enthalpy(c, th, ph, T, Tref, [])
    return tot
