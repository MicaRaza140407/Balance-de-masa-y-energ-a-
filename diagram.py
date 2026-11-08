"""Diagrama de proceso: equipos, corrientes con sus variables y reacciones (matplotlib)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from .render import num

COLORS = {"reactor": "#fde8d0", "mixer": "#dcecf9", "splitter": "#dcecf9", "separator": "#e3f4dc", "generic": "#eeeeee"}


def _label(problem, values, s):
    lines = []
    for c in problem.streams[s]["components"]:
        sym = problem.var[(s, c)]
        tex = problem.reg.tex(sym)
        if values is not None and sym in values:
            lines.append(f"${tex}={num(float(values[sym]), 4)}$ {problem.unit_str}")
        else:
            lines.append(f"${tex}$ [{problem.unit_str}]")
    return lines


def draw(problem, values=None, energy=None):
    units = problem.units
    n = max(len(units), 1)
    fig_w = max(9.0, 4.4 * n + 3.0)
    fig, ax = plt.subplots(figsize=(fig_w, 5.6))
    ax.set_aspect("equal")
    ax.axis("off")
    pos = {}
    bw, bh, gap = 1.9, 1.3, 4.8
    for i, u in enumerate(units):
        x = i * gap
        pos[u["name"]] = (x, 0)
        ax.add_patch(FancyBboxPatch((x - bw / 2, -bh / 2), bw, bh, boxstyle="round,pad=0.05",
                                    fc=COLORS.get(u.get("type", "generic"), "#eeeeee"), ec="#333", lw=1.6))
        ax.text(x, 0.12, u["name"], ha="center", va="center", fontsize=11, weight="bold")
        ax.text(x, -0.22, f"({u.get('type', 'generic')})", ha="center", va="center", fontsize=8, color="#555")
        rxs = [r for r in problem.reactions if r.get("unit") == u["name"]]
        for k, r in enumerate(rxs):
            txt = r["equation"].replace("->", "→")
            if values is not None and problem.xi[r["id"]] in values:
                txt += f"\nξ = {num(float(values[problem.xi[r['id']]]), 4)} {problem.reg.info[problem.xi[r['id']]]['unit']}"
            ax.text(x, -bh / 2 - 0.25 - 0.55 * k, txt, ha="center", va="top", fontsize=8.5, style="italic",
                    bbox=dict(boxstyle="round", fc="#fff8e1", ec="#c9a227", lw=0.8))
        if energy and u["name"] in energy:
            q = energy[u["name"]]["Q"]
            ax.annotate("", xy=(x, bh / 2 + 0.05), xytext=(x, bh / 2 + 0.55) if q >= 0 else (x, bh / 2 + 0.05),
                        arrowprops=dict(arrowstyle="->", color="#c0392b", lw=2))
            ax.text(x + 0.1, bh / 2 + 0.62, f"Q = {num(q, 5)} kJ/{problem.time_u}", ha="left", va="bottom",
                    fontsize=8.5, color="#c0392b")
    src, dst = {}, {}
    for u in units:
        for s in u["outputs"]: src[str(s)] = u["name"]
        for s in u["inputs"]: dst[str(s)] = u["name"]
    idx = {u["name"]: i for i, u in enumerate(units)}
    in_ports, out_ports = {}, {}
    for u in units:
        ins = [str(s) for s in u["inputs"]]
        outs = [str(s) for s in u["outputs"]]
        for k, s in enumerate(ins):
            in_ports[s] = (pos[u["name"]][0] - bw / 2, _spread(k, len(ins), bh))
        for k, s in enumerate(outs):
            out_ports[s] = (pos[u["name"]][0] + bw / 2, _spread(k, len(outs), bh))
    for s in problem.streams:
        arrow = dict(arrowstyle="-|>", mutation_scale=16, lw=1.8, color="#1f4e79")
        lab = "\n".join(_label(problem, values, s))
        if s not in src:                       # alimentación externa
            x1, y1 = in_ports[s]
            a, b = (x1 - 2.2, y1), (x1, y1)
            ax.add_patch(FancyArrowPatch(a, b, **arrow))
            ax.text(a[0] + 1.1, y1 + 0.12, f"({s}) {problem.streams[s].get('name', '')}", ha="center", va="bottom", fontsize=8, weight="bold")
            ax.text(a[0] + 1.1, y1 - 0.12, lab, ha="center", va="top", fontsize=7.5)
        elif s not in dst:                     # producto externo
            x1, y1 = out_ports[s]
            a, b = (x1, y1), (x1 + 2.2, y1)
            ax.add_patch(FancyArrowPatch(a, b, **arrow))
            ax.text(a[0] + 1.1, y1 + 0.12, f"({s}) {problem.streams[s].get('name', '')}", ha="center", va="bottom", fontsize=8, weight="bold")
            ax.text(a[0] + 1.1, y1 - 0.12, lab, ha="center", va="top", fontsize=7.5)
        else:                                  # corriente interna
            a, b = out_ports[s], in_ports[s]
            if idx[dst[s]] > idx[src[s]]:
                rad = 0.0 if idx[dst[s]] - idx[src[s]] == 1 else 0.35
                ax.add_patch(FancyArrowPatch(a, b, connectionstyle=f"arc3,rad={rad}", **arrow))
                ax.text((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + (0.15 if rad == 0 else 0.9), f"({s}) {problem.streams[s].get('name', '')}",
                        ha="center", va="bottom", fontsize=8, weight="bold")
                ax.text((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + (0.0 if rad == 0 else 0.75), lab, ha="center", va="top", fontsize=7.5)
            else:                              # recirculación (hacia atrás): por debajo de los equipos
                k = idx[src[s]] - idx[dst[s]]
                drop = -2.4 - 0.7 * min(k, 3)
                sx = pos[src[s]][0] + bw / 4
                dx = pos[dst[s]][0] - bw / 4
                col = "#7b3f00"
                ax.plot([sx, sx, dx], [-bh / 2 - 0.02, drop, drop], color=col, lw=1.8)
                ax.add_patch(FancyArrowPatch((dx, drop), (dx, -bh / 2 - 0.03), arrowstyle="-|>", mutation_scale=16, lw=1.8, color=col))
                ax.text((sx + dx) / 2, drop - 0.12, f"({s}) {problem.streams[s].get('name', '')} (recirculación)\n" + "   ".join(_label(problem, values, s)),
                        ha="center", va="top", fontsize=7.5, color=col)
    ax.set_title(problem.title, fontsize=11)
    ax.autoscale_view()
    xs = [p[0] for p in pos.values()]
    ax.set_xlim(min(xs) - bw / 2 - 2.5, max(xs) + bw / 2 + 2.7)
    ax.set_ylim(-4.6, 2.2)
    fig.tight_layout()
    return fig


def _spread(k, n, h):
    if n == 1:
        return 0.0
    span = h * 0.7
    return span / 2 - k * span / (n - 1)
