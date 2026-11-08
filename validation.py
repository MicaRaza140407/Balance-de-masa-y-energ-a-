"""Verificación de resultados (paso final de Felder: ¿tiene sentido físico la respuesta?)."""
from . import units as U
from .data import ATOMIC


def validate(problem, values, dof, inconsistent, redundant, energy=None):
    """Devuelve (checks, issues).  checks: lista de (ok, texto). issues: advertencias/errores."""
    checks, issues = [], []
    reg = problem.reg
    # 1) datos faltantes / incoherentes
    free = [s for s in reg.syms if s not in values]
    if free:
        issues.append("DATOS FALTANTES: el sistema queda indeterminado (grados de libertad = "
                      f"{dof['dof']}). Variables sin resolver: {', '.join(reg.plain(s) for s in free)}. "
                      "Agregue datos (flujos, composiciones, conversiones o relaciones) hasta que los grados de libertad sean 0.")
        checks.append((False, "Todas las incógnitas quedaron determinadas"))
    else:
        checks.append((True, "Todas las incógnitas quedaron determinadas (grados de libertad = 0)"))
    if inconsistent:
        for e, x in inconsistent:
            issues.append(f"DATOS INCONSISTENTES: la ecuación {e.id} ({e.why}) queda como 0 = {float(x):.6g} "
                          "tras sustituir; revise los datos del enunciado (p. ej. fracciones que no suman 1).")
        checks.append((False, "Los datos son consistentes entre sí"))
    else:
        checks.append((True, "Los datos son consistentes entre sí" + (f" ({len(redundant)} ecuación(es) redundante(s) satisfecha(s))" if redundant else "")))
    if free:
        return checks, issues
    val = {s: float(v) for s, v in values.items()}
    # 2) residuos de los balances
    worst = 0.0
    for e in problem.eqs:
        r = float(e.expr.subs({s: values[s] for s in e.expr.free_symbols}))
        worst = max(worst, abs(r))
    checks.append((worst < 1e-8, f"Todas las ecuaciones (balances y datos) se cumplen: residuo máximo = {worst:.2e}"))
    # 3) negativos
    neg = [(reg.plain(s), v) for s, v in val.items() if reg.info[s]["kind"] != "extent" and v < -1e-9]
    if neg:
        issues.append("RESULTADOS NO FÍSICOS: flujos negativos → " + "; ".join(f"{n} = {v:.4g}" for n, v in neg)
                      + ". Revise los datos o las condiciones del proceso (p. ej. conversiones o composiciones imposibles).")
    checks.append((not neg, "Ningún flujo de componente es negativo"))
    for r in problem.reactions:
        xi = val[problem.xi[r["id"]]]
        if xi < -1e-9:
            issues.append(f"El avance de la reacción {r['id']} es negativo ({xi:.4g}): la reacción iría en sentido inverso.")
    # 4) fracciones
    bad = []
    for s, info in problem.streams.items():
        tot = sum(val[problem.var[(s, c)]] * U.amount_factor(problem.int_u, "kg", problem.comps[c]["mw"])
                  for c in info["components"])
        if tot > 0:
            sx = sum(val[problem.var[(s, c)]] * U.amount_factor(problem.int_u, "kg", problem.comps[c]["mw"]) / tot
                     for c in info["components"])
            if abs(sx - 1) > 1e-9:
                bad.append(s)
    checks.append((not bad, "Las fracciones de cada corriente suman 1"))
    # 5) balance global de masa y atómicos
    ext_in, ext_out = _external(problem)
    def mass(streams):
        return sum(val[problem.var[(s, c)]] * U.amount_factor(problem.int_u, "kg", problem.comps[c]["mw"])
                   for s in streams for c in problem.streams[s]["components"])
    mi, mo = mass(ext_in), mass(ext_out)
    ok = abs(mi - mo) <= 1e-8 * max(1.0, abs(mi))
    checks.append((ok, f"Balance global de masa: entra {mi:.6g} kg/{problem.time_u}, sale {mo:.6g} kg/{problem.time_u}"))
    if not ok:
        issues.append("El balance global de masa no cierra: revise estequiometría (¿ecuación balanceada? ¿masas molares?).")
    if all(problem.comps[c]["elements"] for c in problem.comps):
        els = set()
        for c in problem.comps.values():
            els |= set(c["elements"])
        okat, txt = True, []
        for el in sorted(els):
            def atoms(streams):
                return sum(val[problem.var[(s, c)]] * U.amount_factor(problem.int_u, "mol", problem.comps[c]["mw"])
                           * problem.comps[c]["elements"].get(el, 0) for s in streams for c in problem.streams[s]["components"])
            a_in, a_out = atoms(ext_in), atoms(ext_out)
            good = abs(a_in - a_out) <= 1e-8 * max(1.0, abs(a_in))
            okat &= good
            txt.append(f"{el}: {a_in:.6g} → {a_out:.6g} mol/{problem.time_u}")
        checks.append((okat, "Balance atómico global (entra → sale): " + "; ".join(txt)))
        if not okat:
            issues.append("El balance atómico no cierra: la reacción no está balanceada o falta/ sobra una especie.")
    # 6) unidades
    checks.append((True, f"Unidades consistentes: todos los datos se convirtieron a {problem.unit_str} antes de operar"))
    # 7) energía
    if energy:
        for u, r in energy.items():
            ok = abs((r["Q"] - r["Ws"]) - r["dH"]) < 1e-8 * max(1.0, abs(r["dH"]))
            checks.append((ok, f"Energía «{u}»: Q − Ws = ΔH → {r['Q'] - r['Ws']:.6g} = {r['dH']:.6g} kJ/{problem.time_u}"))
    return checks, issues


def _external(problem):
    outs_of = {str(s) for u in problem.units for s in u["outputs"]}
    ins_of = {str(s) for u in problem.units for s in u["inputs"]}
    ext_in = [s for s in problem.streams if s not in outs_of]
    ext_out = [s for s in problem.streams if s not in ins_of]
    return ext_in, ext_out
