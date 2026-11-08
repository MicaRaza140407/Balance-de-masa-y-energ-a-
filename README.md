# ⚗️ Balances de Masa y Energía — herramienta educativa (Felder & Rousseau)

Programa en Python que plantea y resuelve **paso a paso** problemas de balance de materia y energía de
Ingeniería Química, con diagrama de proceso, análisis de grados de libertad, desarrollo matemático con unidades
y verificación de resultados. Pensado como herramienta de estudio, no como una caja negra de respuestas.

Metodología basada en:
- Felder & Rousseau, *Principios elementales de los procesos químicos*, 3.ª ed.
- Felder, Rousseau & Bullard, *Elementary Principles of Chemical Processes*, 4th ed.

> Implementa el procedimiento general de esos libros (diagrama rotulado → base de cálculo → grados de libertad →
> balances → resolución → verificación; balances de especies con avance de reacción; Q − Ws = ΔH con estados de
> referencia). **No contiene texto de los libros.** Los datos termodinámicos de `balance_solver/data.py`
> (Cp, ΔH_v, ΔH_f) deben verificarse contra las tablas B.1 y B.2 de su edición; pueden sobrescribirse en el JSON.

## Qué hace

| Requisito | Dónde / cómo |
|---|---|
| Ingreso del enunciado | Pestaña 1 de la interfaz (`app.py`) |
| Identificación de datos, unidades, componentes, condiciones | `parser.analizar_enunciado` (reglas, sin red) y, opcionalmente, interpretación con IA |
| Resolución paso a paso | `mass_balance.eliminate`: cada ecuación, sustitución, despeje y resultado con unidades |
| Diagramas | `diagram.draw`: equipos, corrientes rotuladas con sus variables, reacciones y reciclos |
| Notación | `m` masa, `n` moles, subíndices `m_{i,s}` (componente i, corriente s), `ξ` avance |
| Balances | globales y por componente, con/sin reacción, múltiples unidades, reciclo, derivación, purga, divisores |
| Energía | calor sensible (Cp polinomial), calor latente, cambio de fase, reacción (ΔH°r desde ΔH_f), Q, Ws, T de salida desconocida |
| Validación | datos faltantes (GL>0), datos inconsistentes, balances atómico y de masa globales, flujos negativos, fracciones, unidades, Q − Ws = ΔH |

## Instalación

```bash
git clone <su-repositorio>
cd balance-quimica
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```
Probado con Python 3.12 (debería funcionar en 3.9 o superior).

## Uso

**Interfaz gráfica**
```bash
streamlit run app.py
```
1. Barra lateral → cargue un ejemplo, o pegue su enunciado en la pestaña **Enunciado**.
2. *Identificar datos* muestra lo que se reconoce por reglas.
3. Complete/edite la **especificación JSON** (pestaña 2). Con una `ANTHROPIC_API_KEY` el botón
   *Convertir a especificación con IA* la genera a partir del enunciado; **revísela siempre**.
4. **Solución** → *Resolver*: diagrama, informe completo y descarga en Markdown.

**Línea de comandos**
```bash
python -m balance_solver examples/03_reciclo_reactor.json --out informe.md --png diagrama.png
python -m balance_solver --texto "Se alimentan 1000 kg/h de ..."     # solo identificación de datos
```
**Como librería**
```python
import json
from balance_solver import solve_problem
r = solve_problem(json.load(open("examples/02_combustion_metano.json")))
print(r.value(3, "CO2"), r.energy["Horno"]["Q"], r.ok)
print(r.markdown)
```

## Formato de la especificación (JSON)

```jsonc
{
  "title": "...", "statement": "...", "flow_unit": "kg/h",          // kg, g, lb, ton, mol, kmol, lbmol  /  s, min, h, d
  "components": [{"name": "Etanol", "formula": "C2H5OH", "mw": 46.07 /*opcional*/}],
  "streams": {
    "1": {"name": "Alimentación", "components": ["Etanol","Agua"],  // opcional; por defecto todos
          "total": "1000 kg/h",                                     // flujo total (puede estar en otra base)
          "flows": {"Etanol": "400 kg/h"},                          // flujos de componentes conocidos
          "fractions": {"Etanol": 0.4}, "fraction_basis": "mass"}   // "mass" | "mole"
  },
  "units": [{"name": "Columna", "type": "separator", "inputs": [1], "outputs": [2,3]}],
  "reactions": [{"id":"R1","unit":"Reactor","equation":"A -> B",
                 "conversion":{"reactant":"A","value":0.2,"of_streams":[...] /*opcional*/} /* o */ "extent":"10 mol/h"}],
  "ratios":  [{"a":{"stream":2,"comp":"O2"},"b":{"stream":1,"comp":"CH4"},"value":2.4,"text":"..."}],
  "splits":  [{"unit":"Divisor","stream":5,"fraction":0.1}],        // fracción de la alimentación que va a esa salida (purga/bypass)
  "linear":  [{"terms":[{"stream":3,"comp":"A","coef":1}],"rhs":0,"text":"..."}],
  "basis":   {"stream":1,"total":"100 mol/h"},                      // opcional; si faltan datos de escala se toma 100 en la 1.ª alimentación
  "energy":  {"Tref":25,
              "streams":{"1":{"T":25,"phases":{"CH4":"gas"}}},      // T en °C
              "units":{"Reactor":{"Ws":0}},
              "solve_T":{"unit":"U","stream":2,"Q":25112.4,"bracket":[0,99]},
              "thermo":{"C2H5OH":{"cp":{"liquid":[a,b,c,d]},"ref_phase":"liquid","Hf":-277.6}}}   // opcional: añadir/sobrescribir especies
}
```
Un estado estacionario se asume siempre. Cp [kJ/(mol·°C)] = a + bT + cT² + dT³ (T en °C).
Los ejemplos de `examples/` (destilación, combustión, reciclo con reactor, calor sensible, cambio de fase, T desconocida) sirven de plantilla.

## Estructura

```
app.py                      interfaz Streamlit
balance_solver/
  units.py                  unidades y conversiones masa↔moles↔tiempo
  data.py                   masas atómicas, fórmulas, datos termodinámicos
  mass_balance.py           ecuaciones, grados de libertad, resolución paso a paso
  energy_balance.py         entalpías, calores de reacción, Q − Ws = ΔH
  diagram.py                diagramas de proceso (matplotlib)
  validation.py             verificaciones y detección de errores
  parser.py                 interpretación del enunciado (reglas + IA opcional)
  solver.py                 orquestación e informe Markdown
  render.py                 formato de números y ecuaciones LaTeX
examples/  tests/
```

## Pruebas
```bash
python -m unittest discover -s tests -v
```

## Limitaciones (importante)

- **El lenguaje natural libre no se resuelve solo con reglas.** Las reglas identifican datos, pero armar el flujograma
  (corrientes y unidades) de un enunciado arbitrario requiere interpretación: use la especificación JSON o la
  conversión con IA (que debe revisarse).
- Los sistemas deben ser **lineales** en los flujos por componente: no se resuelven restricciones bilineales
  (p. ej. un divisor con fracción de división desconocida y composición desconocida) ni equilibrios, ecuaciones de estado o cinética.
- Energía: sistemas abiertos en estado estacionario, ΔEc y ΔEp despreciables, sin efecto de presión; fases líquido/gas.
  No incluye sólidos, mezclas no ideales (calores de solución/mezclado), humedad ni vapor saturado por tablas.
- Estado no estacionario y procesos por lotes con energía: no incluidos.

## Cómo ampliar

- **Nuevas especies:** agregue la entrada en `data.THERMO` o en `energy.thermo` del JSON.
- **Nuevo tipo de relación:** agregue un bloque en `MassBalanceProblem._make_equations` que llame a `self._add(id, expr, categoría, tex, explicación)`; el resto (grados de libertad, solución, verificación) funciona sin cambios.
- **Nuevas verificaciones:** `validation.validate`.
