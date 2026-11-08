"""Solucionador educativo de balances de materia y energía (metodología de Felder & Rousseau)."""
from .solver import solve_problem, Result

__all__ = ["solve_problem", "Result"]
__version__ = "1.0.0"
