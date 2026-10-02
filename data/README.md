# Data

`raw/dataset_limpio.xlsx` — 900 records of cyclic triaxial resilient-modulus tests on a granular subbase material
modified with sugarcane bagasse ash (SCBA).

**Experimental structure.** One physical specimen per SCBA content (0 %, 5 %, 10 %) × 15 loading sequences × 20 consecutive
readings per sequence = 900 records. The 20 readings of a sequence are successive readings of the same specimen (a time series),
**not** independent specimens. `Experimento` identifies only the SCBA level.

**Provenance.** Cyclic triaxial resilient-modulus tests, Universidad Cooperativa de Colombia.
Test standard / protocol: [COMPLETA, p. ej. norma y número de secuencias]. Laboratory: [COMPLETA]. Test dates: [COMPLETA]. Instrumentation (LVDTs, load cell): [COMPLETA].

**License.** The data are released under **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/); the code is released under MIT (see `LICENSE`).

**Data-quality screening (notebook, section 4.0).** 25 of the 900 readings are excluded from modelling:
20 readings (10 % SCBA, sequence 15, N = 881–900) with zero axial resilient deformation (sensor fault), and 5 readings
(10 % SCBA, sequence 1, N = 601–605) whose cyclic stress is more than 15 % away from the median of the block (seating cycles).
The raw file is left unchanged.

## Data dictionary

| Column | Unit | Description | Role in the study |
|---|---|---|---|
| `N` | – | Reading index (1–900) | metadata (never a predictor) |
| `Experimento` | – | SCBA content label (0 %, 5 %, 10 %) = specimen | metadata → `% CBCA` |
| `Sequence Number` | – | Loading sequence (1–15) | metadata / validation groups |
| `Axial Resilient Modulus-MPa` | MPa | Resilient modulus from the axial LVDTs, Mr = Δσd / εr | **target** |
| `Actuator Resilient Modulus-MPa` | MPa | Resilient modulus from actuator displacement | excluded (leakage) |
| `Axial Resilient def.-mm`, `Axial1/2 Resilient Def-mm`, `Av Axial Reilient Def-mm` | mm | Recoverable axial deformation (the typo "Reilient" is kept from the source file) | excluded (leakage) |
| `Axial Resilient ustrain-µE`, `Axial1/2 Resilient ustrain-µE`, `Actuator Resilient ustrain-µE` | µε | Recoverable axial strain | excluded (leakage) |
| `Axial Permanent def.-mm`, `Axial1/2 Permanent def.-mm`, `Actuator permanent def.-mm` | mm | Accumulated permanent deformation | M2/M3 only (measured during the test) |
| `Axial Permanent strain-%`, `Actuator permanent strain-%` | % | Accumulated permanent strain | exploratory only |
| `Cyclic Axial Load-kN`, `Contact Load-kN` | kN | Applied loads | not used (redundant with stresses) |
| `Cyclic Axial Stress-kPa` | kPa | Cyclic (deviator) stress, Δσd | predictor |
| `Contact Stress-kPa` | kPa | Contact stress | predictor (M3, M4) |
| `Confining Stress-kPa` | kPa | Confining stress, σ3 | predictor |

Derived in the notebook: `% CBCA` (0, 5, 10), `theta = 3·σ3 + Δσd` (bulk stress) and `Block` (SCBA level × sequence, 45 blocks).
