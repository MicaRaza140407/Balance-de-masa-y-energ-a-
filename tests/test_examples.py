import glob, json, os, unittest, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from balance_solver import solve_problem
from balance_solver.diagram import draw

EX = os.path.join(os.path.dirname(__file__), "..", "examples")


def load(n):
    with open(os.path.join(EX, n), encoding="utf-8") as f:
        return json.load(f)


class TestExamples(unittest.TestCase):
    def test_all_examples_verify(self):
        for f in sorted(glob.glob(os.path.join(EX, "*.json"))):
            with open(f, encoding="utf-8") as fh:
                r = solve_problem(json.load(fh))
            self.assertTrue(r.ok, f"{os.path.basename(f)}: {r.issues} {[c for c in r.checks if not c[0]]}")
            draw(r.problem, r.values, r.energy)

    def test_distillation(self):
        r = solve_problem(load("01_destilacion.json"))
        d = r.value(2, "Etanol") + r.value(2, "Agua")
        b = r.value(3, "Etanol") + r.value(3, "Agua")
        self.assertAlmostEqual(d, 350 / 0.85, places=6)
        self.assertAlmostEqual(b, 1000 - 350 / 0.85, places=6)

    def test_combustion(self):
        r = solve_problem(load("02_combustion_metano.json"))
        self.assertAlmostEqual(r.value(3, "CO2"), 100)
        self.assertAlmostEqual(r.value(3, "H2O"), 200)
        self.assertAlmostEqual(r.value(3, "O2"), 40)
        self.assertAlmostEqual(r.value(2, "N2"), 240 * 79 / 21, places=6)
        self.assertAlmostEqual(r.energy["Horno"]["react"], -890.33 * 100, places=3)
        self.assertLess(r.energy["Horno"]["Q"], 0)

    def test_recycle(self):
        r = solve_problem(load("03_reciclo_reactor.json"))
        self.assertAlmostEqual(r.value(2, "A"), 500)
        self.assertAlmostEqual(r.value(5, "A"), 400)
        self.assertAlmostEqual(r.value(4, "B"), 100)

    def test_sensible_heat(self):
        r = solve_problem(load("04_calentamiento_agua.json"))
        n = 100000 / 18.015
        self.assertAlmostEqual(r.energy["Calentador"]["Q"], n * 0.0754 * 60, delta=1e-6 * n)

    def test_phase_change(self):
        r = solve_problem(load("05_agua_a_vapor.json"))
        n = 36000 / 18.015
        q = r.energy["Caldera"]["Q"] / n
        self.assertTrue(47.9 < q < 48.2, q)   # 5.655 + 40.656 + ~1.72 kJ/mol

    def test_solve_temperature(self):
        r = solve_problem(load("06_temperatura_de_salida.json"))
        self.assertIn("T_{2} = 80", r.markdown)

    def test_missing_data_detected(self):
        spec = load("01_destilacion.json")
        del spec["streams"]["3"]["fractions"]
        r = solve_problem(spec)
        self.assertFalse(r.ok)
        self.assertTrue(any("FALTANTES" in i for i in r.issues))

    def test_inconsistent_data_detected(self):
        spec = load("01_destilacion.json")
        spec["streams"]["1"]["fractions"] = {"Etanol": 0.4, "Agua": 0.7}
        r = solve_problem(spec)
        self.assertTrue(any("INCONSISTENTES" in i for i in r.issues))

    def test_unbalanced_reaction_detected(self):
        spec = load("02_combustion_metano.json")
        spec["reactions"][0]["equation"] = "CH4 + O2 -> CO2 + 2 H2O"
        r = solve_problem(spec)
        self.assertTrue(any("atómico" in i or "masa" in i for i in r.issues))

    def test_auto_basis(self):
        spec = load("01_destilacion.json")
        del spec["streams"]["1"]["total"]
        r = solve_problem(spec)
        self.assertTrue(r.ok)
        self.assertAlmostEqual(r.value(1, "Etanol") + r.value(1, "Agua"), 100)


if __name__ == "__main__":
    unittest.main()
