"""Orquestador: especificación del problema → informe paso a paso, diagrama y verificación."""
from dataclasses import dataclass, field
import copy

from . import units as U
from .mass_balance import MassBalanceProblem, eliminate
from .energy_balance import solve_energy
from .validation import validate
from .render import num


@dataclass
class Result:
    problem: MassBalanceProblem
    values: dict
    markdown: str
    checks: list
    issues: list
    dof: dict
    energy: dict = field(default_factory=dict)

    @property
    def ok(self):
        return all(c[0] for c in self.checks) and not self.issues

    def value(self, stream, comp):
        return float(self.values[self.problem.var[(str(stream), comp)]])

    def extent(self, rid):
        return float(self.values[self.problem.xi[rid]])


def _has_recycle(problem):
    graph = {u["name"]: set() for u in problem.units}
    for u in problem.units:
        for o in u["outputs"]:
            for v in problem.units:
                if str(o) in map(str, v["inputs"]):
                    graph[u["name"]].add(v["name"])
    state = {}
    def dfs(n):
        state[n] = 1
        for m in graph[n]:
            if state.get(m) == 1 or (m not in state and dfs(m)):
                return True
        state[n] = 2
        return False
    return any(dfs(n) for n in graph if n not in state)


METHOD_TEXT = {
    "mixer": "Mezclador: varias entradas, una salida. Se escriben balances de cada componente (entrada = salida).",
    "splitter": "Divisor: las corrientes de salida tienen la misma composición que la alimentación; se usa la fracción de división.",
    "separator": "Separador: una alimentación que se reparte en productos de distinta composición; un balance por componente.",
    "reactor": "Reactor: balance de cada especie molecular con el término ν·ξ (método del avance de reacción).",
    "generic": "Unidad sin reacción: balance de cada especie (entrada = salida).",
}


def build_report(problem, dof, auto_basis_stream, log, values, checks, issues, energy_md):
    reg = problem.reg
    L = []
    L.append(f"# {problem.title}\n")
    if problem.statement:
        L.append(f"> {problem.statement}\n")
    # 1. interpretación
    L.append("## 1. Interpretación del problema\n")
    L.append("**Componentes y masas molares** (calculadas a partir de las masas atómicas):\n")
    L.append("| Componente | Fórmula | Cálculo de M |\n|---|---|---|")
    for n, c in problem.comps.items():
        L.append(f"| {n} | {c['formula']} | {c['mw_steps']} |")
    L.append("")
    L.append(f"**Unidades:** los flujos se expresan en `{problem.spec.get('flow_unit', 'kg/h')}`; las variables internas "
             f"están en `{problem.unit_str}`. {problem.basis_reason}\n")
    L.append("**Notación:** $m_{i,s}$ = flujo de masa del componente $i$ en la corriente $s$; "
             "$n_{i,s}$ = flujo molar; $\\xi_r$ = avance de la reacción $r$.\n")
    kinds = []
    if problem.reactions: kinds.append("con reacción química")
    else: kinds.append("sin reacción química")
    if len(problem.units) > 1: kinds.append("de unidades múltiples")
    if _has_recycle(problem): kinds.append("con recirculación")
    if problem.spec.get("energy"): kinds.append("con balance de energía")
    L.append(f"**Tipo de problema:** {', '.join(kinds)}.\n")
    # 2. metodología
    L.append("## 2. Metodología (Felder & Rousseau)\n")
    L.append("Se sigue el procedimiento general de los capítulos 4 (balances de materia) y 7–9 (energía): "
             "(1) dibujar y rotular el diagrama de flujo, (2) elegir una base de cálculo, (3) hacer el análisis de grados de "
             "libertad, (4) escribir los balances y las ecuaciones de los datos y resolverlos, (5) verificar el resultado.\n")
    seen = set()
    for u in problem.units:
        t = u.get("type", "generic")
        if (t, bool(problem.reactions)) not in seen:
            seen.add((t, bool(problem.reactions)))
            L.append(f"- **{u['name']}** — {METHOD_TEXT.get(t, METHOD_TEXT['generic'])}")
    if problem.reactions:
        L.append("- **Por qué el avance de reacción:** con una o varias reacciones conocidas, ξ permite escribir un solo "
                 "balance por especie, con coeficientes estequiométricos, y relaciona directamente la conversión con ξ.")
    if _has_recycle(problem):
        L.append("- **Recirculación:** el sistema se plantea completo (todas las unidades a la vez); el balance global "
                 "(alrededor de todo el proceso) deja fuera las corrientes internas y sirve como verificación.")
    L.append("")
    L.append("> Nota: este programa implementa el procedimiento general de los libros; las secciones exactas pueden variar según la edición.\n")
    # 3. diagrama
    L.append("## 3. Diagrama de proceso\n")
    L.append("Ver el diagrama generado (corrientes rotuladas con sus variables; los valores aparecen tras resolver).\n")
    # 4. datos
    L.append("## 4. Datos del enunciado por corriente\n")
    for s, info in problem.streams.items():
        bits = []
        if "total" in info: bits.append(f"flujo total = {info['total']}")
        for c, q in info.get("flows", {}).items(): bits.append(f"{c}: {q}")
        for c, x in info.get("fractions", {}).items(): bits.append(f"fracción de {c} = {x}")
        L.append(f"- **Corriente {s}** ({info.get('name', '')}): componentes {', '.join(info['components'])}; "
                 f"{'; '.join(bits) if bits else 'sin datos directos (incógnita)'}")
    L.append("")
    if auto_basis_stream:
        L.append(f"**Base de cálculo:** el enunciado no fija ningún flujo absoluto; se toma 100 {problem.unit_str} en la "
                 f"corriente {auto_basis_stream}.\n")
    elif any("total" in i or i.get("flows") for i in problem.streams.values()) or problem.spec.get("basis"):
        L.append("**Base de cálculo:** se usa como base el flujo dado en el enunciado.\n")
    # 5. variables y GL
    L.append("## 5. Variables incógnita y análisis de grados de libertad\n")
    L.append("| Variable | Significado | Unidad |\n|---|---|---|")
    for s in reg.syms:
        i = reg.info[s]
        if i["kind"] == "extent":
            what = "avance de la reacción"
        else:
            what = f"{'masa' if i['kind']=='mass' else 'moles'} de {i['comp']} en la corriente {i['stream']}"
        L.append(f"| ${i['tex']}$ | {what} | {i['unit']} |")
    L.append("")
    L.append(f"Variables desconocidas: **{dof['unknowns']}**. Ecuaciones planteadas: **{dof['equations']}** "
             f"(independientes: **{dof['rank']}**; redundantes: {dof['redundant']}, p. ej. fracciones que suman 1).\n")
    L.append(f"$$\\text{{GL}} = N_{{\\text{{incógnitas}}}} - N_{{\\text{{ecuaciones independientes}}}} = {dof['unknowns']} - {dof['rank']} = {dof['dof']}$$")
    L.append("Los datos específicos del enunciado se cuentan como ecuaciones (formulación equivalente a contar solo las incógnitas "
             "rotuladas). " + ("GL = 0: el problema está **especificado correctamente** y tiene solución única.\n"
                               if dof["dof"] == 0 else
                               ("GL > 0: faltan datos (problema subespecificado).\n" if dof["dof"] > 0 else "GL < 0: hay datos de más.\n")))
    # 6. ecuaciones
    L.append("## 6. Ecuaciones\n")
    for e in problem.eqs:
        L.append(f"**{e.id}** ({e.category}) — {e.why}\n\n$$ {e.tex} $$\n")
    # 7. resolución
    L.append("## 7. Resolución paso a paso\n")
    L.extend(x + "\n" for x in log)
    # 8. resultados
    L.append("## 8. Resultados\n")
    if all(s in values for s in reg.syms):
        L.append(results_table(problem, values))
        for r in problem.reactions:
            L.append(f"\nAvance de reacción: $\\xi_{{{r['id']}}} = {num(float(values[problem.xi[r['id']]]))}$ {reg.info[problem.xi[r['id']]]['unit']}")
        L.append("")
    else:
        L.append("El problema no pudo resolverse completamente (ver verificación).\n")
    if energy_md:
        L.extend(energy_md)
    # 9. verificación
    L.append("\n## 9. Verificación de resultados\n")
    for ok, t in checks:
        L.append(f"- {'✅' if ok else '❌'} {t}")
    L.append("")
    if issues:
        L.append("**Advertencias / errores detectados:**\n")
        for i in issues: L.append(f"- ⚠️ {i}")
    else:
        L.append("No se detectaron errores ni inconsistencias.")
    return "\n".join(L)


def results_table(problem, values):
    tu = problem.time_u
    rows = ["| Corriente | Componente | " + ("n" if problem.kind == "mole" else "m") + f" [{problem.unit_str}] | "
            f"m [kg/{tu}] | fracción másica | fracción molar |", "|---|---|---|---|---|---|"]
    for s, info in problem.streams.items():
        data = []
        for c in info["components"]:
            v = float(values[problem.var[(s, c)]])
            mw = problem.comps[c]["mw"]
            data.append((c, v, v * U.amount_factor(problem.int_u, "kg", mw), v * U.amount_factor(problem.int_u, "mol", mw)))
        tm, tmol, tv = sum(d[2] for d in data), sum(d[3] for d in data), sum(d[1] for d in data)
        for c, v, m, mol in data:
            rows.append(f"| {s} | {c} | {num(v)} | {num(m)} | {num(m / tm) if tm else '0'} | {num(mol / tmol) if tmol else '0'} |")
        rows.append(f"| **{s}** | **total** | **{num(tv)}** | **{num(tm)}** | 1 | 1 |")
    return "\n".join(rows)


def solve_problem(spec):
    """Punto de entrada: recibe el diccionario del problema y devuelve un `Result`."""
    spec = copy.deepcopy(spec)
    problem = MassBalanceProblem(spec)
    dof = problem.dof()
    auto = None
    if not dof["has_scale"] and spec.get("auto_basis", True):
        auto = problem.add_auto_basis(100)
        dof = problem.dof()
    log = []
    values, inconsistent, redundant = eliminate(problem, log)
    energy, energy_md, issues_e = {}, [], []
    if spec.get("energy") and all(s in values for s in problem.reg.syms) and not inconsistent:
        try:
            energy = solve_energy(problem, values, energy_md, issues_e)
        except ValueError as e:
            issues_e.append(f"Balance de energía: {e}")
    elif spec.get("energy"):
        issues_e.append("No se pudo hacer el balance de energía porque el balance de materia no quedó resuelto.")
    checks, issues = validate(problem, values, dof, inconsistent, redundant, energy)
    issues += issues_e
    md = build_report(problem, dof, auto, log, values, checks, issues, energy_md)
    return Result(problem, values, md, checks, issues, dof, energy)
