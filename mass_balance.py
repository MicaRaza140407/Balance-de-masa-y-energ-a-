"""Balances de materia (Felder & Rousseau, caps. 4): formulación y resolución paso a paso.

Procedimiento (Felder): 1) diagrama de flujo con variables rotuladas, 2) base de cálculo,
3) análisis de grados de libertad, 4) balances + ecuaciones de datos, resolución.
Todas las ecuaciones del problema son LINEALES en las cantidades por componente
(n_{i,s} o m_{i,s}) y en los avances de reacción (ξ), por lo que se resuelve por
despeje y sustitución mostrando cada operación.
"""
import re
from dataclasses import dataclass
import sympy as sp

from . import units as U
from .data import molar_mass, parse_formula
from .render import num, lin_tex

R = sp.Rational


def rat(x):
    return R(str(float(x))) if not isinstance(x, sp.Basic) else x


@dataclass
class Eq:
    id: str
    expr: sp.Expr
    category: str      # dato | balance | reaccion | relacion | base
    tex: str           # forma natural de la ecuación
    why: str           # explicación (de dónde sale y por qué se aplica)


class Registry:
    """Relaciona cada símbolo con su nombre, LaTeX y unidad."""

    def __init__(self):
        self.syms, self.info = [], {}

    def add(self, name, tex, plain, unit, kind, comp=None, stream=None):
        s = sp.Symbol(name)
        self.syms.append(s)
        self.info[s] = dict(tex=tex, plain=plain, unit=unit, kind=kind, comp=comp, stream=stream)
        return s

    def tex(self, s): return self.info[s]["tex"]
    def plain(self, s): return self.info[s]["plain"]
    def order(self, s): return self.syms.index(s)

    def order_of(self, term):
        fs = sorted(term.free_symbols, key=self.order)
        return self.order(fs[0]) if fs else 10 ** 6


class MassBalanceProblem:
    def __init__(self, spec):
        self.spec = spec
        self.title = spec.get("title", "Problema de balance de materia")
        self.statement = spec.get("statement", "")
        self.amount_u, self.time_u = U.parse_flow_unit(spec.get("flow_unit", "kg/h"))
        self.comps = {}
        self._read_components()
        self.reactions = self._read_reactions()
        if self.reactions:
            self.int_u = U.MOLE_OF[self.amount_u]
            self.basis_reason = ("Hay reacción química: Felder trabaja las reacciones en términos de moles "
                                 "(coeficientes estequiométricos y avance de reacción ξ), por lo que las variables "
                                 "se definen en moles (n) y los datos en masa se convierten con las masas molares.")
        else:
            self.int_u = self.amount_u
            self.basis_reason = ("No hay reacción química: los balances de especies moleculares se plantean "
                                 "directamente en la base en que están los datos.")
        self.kind = U.kind_of(self.int_u)
        self.p = "m" if self.kind == "mass" else "n"
        self.unit_str = f"{self.int_u}/{self.time_u}"
        self.streams = {str(k): v for k, v in spec["streams"].items()}
        self.units = spec.get("units", [])
        self.reg = Registry()
        self.var = {}          # (stream, comp) -> symbol
        self.xi = {}           # reaction id -> symbol
        self._make_vars()
        self.eqs = []
        self._make_equations()

    # ---------------------------------------------------------- lectura
    def _read_components(self):
        for c in self.spec["components"]:
            d = {"name": c} if isinstance(c, str) else dict(c)
            name = d["name"]
            d.setdefault("formula", name)
            d["elements"] = None
            try:
                d["elements"] = parse_formula(d["formula"])
            except ValueError:
                pass
            if "mw" in d:
                d["mw_steps"] = f"M = {d['mw']} g/mol (dato del enunciado)"
            else:
                try:
                    d["mw"], d["mw_steps"] = molar_mass(d["formula"])
                except ValueError as e:
                    raise ValueError(f"Componente '{name}': {e}. Indique su masa molar con \"mw\".")
            self.comps[name] = d

    def _read_reactions(self):
        out = []
        for r in self.spec.get("reactions", []):
            st = {}
            lhs, rhs = re.split(r"->|→|=>", r["equation"])
            for side, sgn in ((lhs, -1), (rhs, 1)):
                for term in side.split("+"):
                    m = re.match(r"^\s*(\d*\.?\d*(?:/\d+)?)\s*(\S.*?)\s*$", term)
                    if not m:
                        raise ValueError(f"Reacción mal escrita: {r['equation']}")
                    coef = m.group(1)
                    if "/" in coef:
                        a, b = coef.split("/"); coef = R(int(a), int(b))
                    else:
                        coef = R(coef) if coef else R(1)
                    sp_name = m.group(2)
                    if sp_name not in self.comps:
                        raise ValueError(f"La especie '{sp_name}' de la reacción no está en 'components'.")
                    st[sp_name] = st.get(sp_name, 0) + sgn * coef
            rr = dict(r)
            rr["stoich"] = st
            out.append(rr)
        return out

    def _make_vars(self):
        for s, info in self.streams.items():
            comps = info.get("components", list(self.comps))
            info["components"] = comps
            for c in comps:
                if c not in self.comps:
                    raise ValueError(f"Corriente {s}: componente '{c}' no declarado en 'components'.")
                tex = f"{self.p}_{{\\mathrm{{{c}}},{s}}}"
                self.var[(s, c)] = self.reg.add(f"{self.p}[{c},{s}]", tex, f"{self.p}_{c},{s}",
                                                self.unit_str, self.kind, c, s)
        for r in self.reactions:
            self.xi[r["id"]] = self.reg.add(f"xi[{r['id']}]", f"\\xi_{{{r['id']}}}", f"xi_{r['id']}",
                                            f"{U.MOLE_OF[self.amount_u] if self.kind=='mass' else self.int_u}/{self.time_u}",
                                            "extent")

    # ---------------------------------------------------------- ayudas
    def total_expr(self, s, w_unit=None):
        """Σ_i w_i·a_{i,s}, con w_i que convierte la unidad interna a w_unit."""
        e = 0
        for c in self.streams[s]["components"]:
            w = U.amount_factor(self.int_u, w_unit or self.int_u, self.comps[c]["mw"]) if w_unit else 1.0
            e += rat(w) * self.var[(s, c)]
        return e

    def total_tex(self, s):
        return " + ".join(self.reg.tex(self.var[(s, c)]) for c in self.streams[s]["components"])

    def term(self, s, comp):
        s = str(s)
        if s not in self.streams:
            raise ValueError(f"Corriente '{s}' no existe.")
        if comp in (None, "__total__"):
            return self.total_expr(s), f"{self.p}_{{{s}}}"
        if (s, comp) not in self.var:
            raise ValueError(f"La corriente {s} no contiene '{comp}'.")
        return self.var[(s, comp)], self.reg.tex(self.var[(s, comp)])

    def _conv_value(self, q, comp_mw=None):
        """Convierte (valor, unidades) a la unidad interna por tiempo del problema."""
        v, au, tu = U.parse_quantity(q)
        au = au or self.amount_u
        tu = tu or self.time_u
        f = U.amount_factor(au, self.int_u, comp_mw) * U.time_factor(tu, self.time_u)
        return v * f, v, au, tu

    def _add(self, id_, expr, cat, tex, why):
        self.eqs.append(Eq(id_, sp.expand(expr), cat, tex, why))

    # ---------------------------------------------------------- ecuaciones
    def _make_equations(self):
        n = 0

        def nid(prefix):
            nonlocal n
            n += 1
            return f"{prefix}{n}"

        # 1) datos de cada corriente
        for s, info in self.streams.items():
            if "total" in info:
                q = U.parse_quantity(info["total"])
                au = q[1] or self.amount_u
                tu = q[2] or self.time_u
                value = q[0] * U.time_factor(tu, self.time_u)
                expr = self.total_expr(s, au) - rat(value)
                self._add(nid("D"), expr, "dato",
                          f"{self.total_tex(s)} = {num(value)}\\ \\text{{{au}/{self.time_u}}}"
                          if U.kind_of(au) == self.kind else
                          f"\\sum_i w_i\\,{self.p}_{{i,{s}}} = {num(value)}\\ \\text{{{au}/{self.time_u}}}",
                          f"Flujo total conocido de la corriente {s} ({q[0]:g} {au}/{tu}); "
                          f"la suma de las cantidades de sus componentes debe igualarlo.")
            for c, qv in info.get("flows", {}).items():
                if (s, c) not in self.var:
                    raise ValueError(f"Corriente {s}: '{c}' no es componente de la corriente.")
                value, v0, au, tu = self._conv_value(qv, self.comps[c]["mw"])
                self._add(nid("D"), self.var[(s, c)] - rat(value), "dato",
                          f"{self.reg.tex(self.var[(s, c)])} = {num(value)}\\ \\text{{{self.unit_str}}}",
                          f"Flujo conocido de {c} en la corriente {s} ({v0:g} {au}/{tu}).")
            fr = info.get("fractions", {})
            if fr:
                fb = info.get("fraction_basis", "mass" if U.kind_of(self.amount_u) == "mass" else "mole")
                fu = "kg" if fb == "mass" else "mol"
                ftxt = "másica" if fb == "mass" else "molar"
                for c, x in fr.items():
                    x = float(str(x).rstrip("%")) / 100 if isinstance(x, str) and x.endswith("%") else float(x)
                    if (s, c) not in self.var:
                        if x != 0:
                            raise ValueError(f"Corriente {s}: se da fracción de '{c}' pero no es componente.")
                        continue
                    tot = sum((rat(U.amount_factor(self.int_u, fu, self.comps[j]["mw"]))
                               * self.var[(s, j)] for j in info["components"]))
                    wi = rat(U.amount_factor(self.int_u, fu, self.comps[c]["mw"]))
                    sym = "x" if fb == "mass" else "y"
                    word = "masa" if fb == "mass" else "moles"
                    self._add(nid("D"), wi * self.var[(s, c)] - rat(x) * tot, "dato",
                              sym + "_{\\mathrm{" + c + "}," + s + "} = \\dfrac{\\text{" + word + " de " + c + "}}{\\text{"
                              + word + " total}} = " + num(x),
                              f"Fracción {ftxt} de {c} en la corriente {s} ({x:g}): "
                              f"cantidad de {c} = fracción × cantidad total de la corriente.")
        # base explícita
        b = self.spec.get("basis")
        if b:
            s = str(b["stream"])
            q = U.parse_quantity(b["total"])
            au, tu = q[1] or self.amount_u, q[2] or self.time_u
            value = q[0] * U.time_factor(tu, self.time_u)
            self._add(nid("B"), self.total_expr(s, au) - rat(value), "base",
                      f"{self.total_tex(s)} = {num(value)}\\ \\text{{{au}/{self.time_u}}}",
                      f"Base de cálculo elegida: {q[0]:g} {au}/{tu} en la corriente {s} (paso 2 de Felder).")

        # 2) balances por unidad
        by_rx = {}
        for r in self.reactions:
            by_rx.setdefault(r.get("unit"), []).append(r)
        self.unit_streams = {}
        for u in self.units:
            ins, outs = [str(x) for x in u["inputs"]], [str(x) for x in u["outputs"]]
            for s in ins + outs:
                if s not in self.streams:
                    raise ValueError(f"La unidad '{u['name']}' usa la corriente {s}, que no está definida.")
            self.unit_streams[u["name"]] = (ins, outs)
            rxs = by_rx.get(u["name"], [])
            comps_u = sorted({c for s in ins + outs for c in self.streams[s]["components"]},
                             key=list(self.comps).index)
            for c in comps_u:
                out_e = sum((self.var[(s, c)] for s in outs if (s, c) in self.var), sp.Integer(0))
                in_e = sum((self.var[(s, c)] for s in ins if (s, c) in self.var), sp.Integer(0))
                gen = sum((self.xi_stoich(r, c) * self.xi[r["id"]] for r in rxs), sp.Integer(0))
                tex_out = " + ".join(self.reg.tex(self.var[(s, c)]) for s in outs if (s, c) in self.var) or "0"
                tex_in = " + ".join(self.reg.tex(self.var[(s, c)]) for s in ins if (s, c) in self.var) or "0"
                tex = f"{tex_out} = {tex_in}"
                why = f"Balance de {c} en «{u['name']}» (estado estacionario): salida = entrada"
                if rxs:
                    terms = " ".join(f"{'+' if self.xi_stoich(r, c) > 0 else '-'} {abs(self.xi_stoich(r, c))}\\,\\xi_{{{r['id']}}}"
                                     for r in rxs if self.xi_stoich(r, c) != 0)
                    tex += f" {terms}"
                    why += " + generación − consumo por reacción (ν·ξ, balance por avance de reacción)"
                why += "." if not rxs else "."
                self._add(nid("E"), out_e - in_e - gen, "balance", tex, why)
            # divisores con fracción de división
            for sp_ in self.spec.get("splits", []):
                if sp_["unit"] != u["name"]:
                    continue
                s_out = str(sp_["stream"])
                f = float(sp_["fraction"])
                if len(ins) != 1:
                    raise ValueError("Un divisor debe tener una sola corriente de entrada.")
                for c in self.streams[s_out]["components"]:
                    if (ins[0], c) not in self.var:
                        continue
                    self._add(nid("S"), self.var[(s_out, c)] - rat(f) * self.var[(ins[0], c)], "relacion",
                              f"{self.reg.tex(self.var[(s_out, c)])} = {num(f)}\\,{self.reg.tex(self.var[(ins[0], c)])}",
                              f"Las corrientes de un divisor tienen la misma composición que la alimentación; la "
                              f"corriente {s_out} recibe la fracción {f:g} de cada componente.")
        # 3) reacciones: conversión / avance
        for r in self.reactions:
            xi = self.xi[r["id"]]
            unit = next((u for u in self.units if u["name"] == r.get("unit")), None)
            if "extent" in r:
                value, v0, au, tu = self._conv_value(r["extent"], None)
                self._add(nid("R"), xi - rat(value), "reaccion", f"\\xi_{{{r['id']}}} = {num(value)}",
                          f"Avance de la reacción {r['id']} dado en el enunciado.")
            conv = r.get("conversion")
            if conv:
                a = conv["reactant"]
                x = float(conv["value"])
                of = [str(z) for z in conv.get("of_streams", unit["inputs"] if unit else [])]
                feed = sum((self.var[(s, a)] for s in of if (s, a) in self.var), sp.Integer(0))
                nu = self.xi_stoich(r, a)
                if nu >= 0:
                    raise ValueError(f"'{a}' no es reactivo en {r['id']}.")
                what = "de paso simple (sobre la alimentación al reactor)" if not conv.get("of_streams") \
                    else f"sobre las corrientes {', '.join(of)}"
                self._add(nid("R"), nu * xi + rat(x) * feed, "reaccion",
                          f"{abs(nu)}\\,\\xi_{{{r['id']}}} = {num(x)}\\,\\left({' + '.join(self.reg.tex(self.var[(s, a)]) for s in of if (s, a) in self.var)}\\right)",
                          f"Conversión {what} de {a} = {x:g}: reactivo consumido (|ν|·ξ) = X × reactivo alimentado.")
        # 4) relaciones adicionales
        for rt in self.spec.get("ratios", []):
            ea, ta = self.term(rt["a"]["stream"], rt["a"].get("comp"))
            eb, tb = self.term(rt["b"]["stream"], rt["b"].get("comp"))
            val = float(rt["value"])
            self._add(nid("Q"), ea - rat(val) * eb, "relacion", f"{ta} = {num(val)}\\,{tb}",
                      rt.get("text", "Relación entre flujos dada en el enunciado."))
        for ln in self.spec.get("linear", []):
            e, tx = 0, []
            for t in ln["terms"]:
                ex, tt = self.term(t["stream"], t.get("comp"))
                e += rat(t.get("coef", 1)) * ex
                tx.append(f"{num(t.get('coef', 1))}\\,{tt}")
            val = float(ln.get("rhs", 0))
            self._add(nid("Q"), e - rat(val), "relacion", f"{' + '.join(tx)} = {num(val)}",
                      ln.get("text", "Relación lineal dada en el enunciado."))

    @staticmethod
    def xi_stoich(r, c):
        return r["stoich"].get(c, sp.Integer(0))

    # ---------------------------------------------------------- grados de libertad
    def dof(self):
        syms = self.reg.syms
        A, b = sp.linear_eq_to_matrix([e.expr for e in self.eqs], syms)
        rk = A.rank()
        Ab = A.row_join(b)
        rk_aug = Ab.rank()
        return dict(unknowns=len(syms), equations=len(self.eqs), rank=rk, dof=len(syms) - rk,
                    redundant=len(self.eqs) - rk, consistent=(rk == rk_aug),
                    has_scale=any(sp.expand(e.expr).subs({s: 0 for s in syms}) != 0 for e in self.eqs))

    def add_auto_basis(self, amount=100):
        feeds = [s for s in self.streams if not any(s in u["outputs"] or str(s) in map(str, u["outputs"])
                                                    for u in self.units)]
        s = feeds[0]
        value = amount
        self._add("B1", self.total_expr(s) - rat(value), "base",
                  f"{self.total_tex(s)} = {num(value)}\\ \\text{{{self.unit_str}}}",
                  f"Los datos no fijan la escala del proceso: se elige como base de cálculo {value} "
                  f"{self.unit_str} en la corriente de alimentación {s} (paso 2 de Felder). "
                  f"Los resultados pueden escalarse proporcionalmente al flujo real.")
        return s


# ---------------------------------------------------------------- resolución
def eliminate(problem, log):
    """Resuelve el sistema por despeje y sustitución, escribiendo cada operación en `log`."""
    reg = problem.reg
    unit = problem.unit_str
    sol, defs = {}, []
    pending = [(e, e.expr) for e in problem.eqs]
    inconsistent, redundant = [], []
    step = 0
    while True:
        nxt = []
        for e, x in pending:
            x = sp.expand(x.subs(sol))
            if not x.free_symbols:
                (inconsistent if abs(float(x)) > 1e-9 else redundant).append((e, x))
            else:
                nxt.append((e, x))
        pending = nxt
        if not pending:
            break
        pending.sort(key=lambda t: (len(t[1].free_symbols), problem.eqs.index(t[0])))
        e, x = pending.pop(0)
        fs = sorted(x.free_symbols, key=reg.order)
        step += 1
        log.append(f"**Paso {step}.** Ecuación **{e.id}** — {e.why}")
        log.append(f"$${e.tex}$$")
        shown = lin_tex(e.expr, reg, known=sol, unit=unit)
        log.append(f"Sustituyendo los valores ya conocidos: $${shown}$$")
        if len(fs) == 1:
            v = fs[0]
            a, c = x.coeff(v), x - x.coeff(v) * v
            val = -c / a
            sol[v] = val
            log.append(f"Queda una sola incógnita. Se despeja: $${reg.tex(v)} = -\\dfrac{{({num(c)})}}{{({num(a)})}} "
                       f"= {num(val)}\\ \\text{{{reg.info[v]['unit']}}}$$")
        else:
            v = fs[-1]
            a, rest = x.coeff(v), x - x.coeff(v) * v
            expr = sp.expand(-rest / a)
            defs.append((v, expr, e))
            pending = [(e2, sp.expand(x2.subs(v, expr))) for e2, x2 in pending]
            log.append(f"Hay {len(fs)} incógnitas: se despeja {reg.tex(v)} en función de las demás y se sustituye "
                       f"en las ecuaciones restantes: $${lin_tex(expr, reg, rhs='', unit=unit).replace(' = ', '')}"
                       f"$$ → $${reg.tex(v)} = {lin_tex(expr, reg, unit=unit, rhs='').replace(' = ', '')}$$")
    for v, expr, e in reversed(defs):
        val = sp.expand(expr.subs(sol))
        if not val.free_symbols:
            sol[v] = val
            step += 1
            log.append(f"**Paso {step} (sustitución regresiva).** Con los valores obtenidos: "
                       f"$${reg.tex(v)} = {lin_tex(expr, reg, known=sol, unit=unit, rhs='').replace(' = ', '')}"
                       f" = {num(val)}\\ \\text{{{reg.info[v]['unit']}}}$$")
    return sol, inconsistent, redundant
