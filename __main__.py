"""Uso:  python -m balance_solver examples/01_destilacion.json [--out informe.md] [--png diagrama.png]
        python -m balance_solver --texto "enunciado ..."   (solo interpretación por reglas)"""
import argparse, json, sys
from . import solve_problem
from .parser import analizar_enunciado


def main(argv=None):
    ap = argparse.ArgumentParser(description="Solucionador de balances de materia y energía (Felder & Rousseau)")
    ap.add_argument("archivo", nargs="?", help="JSON con la especificación del problema")
    ap.add_argument("--out", help="guardar el informe en Markdown")
    ap.add_argument("--png", help="guardar el diagrama de proceso")
    ap.add_argument("--texto", help="analizar un enunciado (extracción por reglas)")
    a = ap.parse_args(argv)
    if a.texto:
        print(analizar_enunciado(a.texto)["markdown"]); return 0
    if not a.archivo:
        ap.error("indique un archivo JSON o --texto")
    with open(a.archivo, encoding="utf-8") as f:
        spec = json.load(f)
    r = solve_problem(spec)
    if a.out:
        open(a.out, "w", encoding="utf-8").write(r.markdown)
    else:
        print(r.markdown)
    if a.png:
        from .diagram import draw
        draw(r.problem, r.values, r.energy).savefig(a.png, dpi=150)
    return 0 if r.ok else 1


if __name__ == "__main__":
    sys.exit(main())
