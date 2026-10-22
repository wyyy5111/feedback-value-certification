# Source provenance

This project packages independently implemented rational value-evaluation algorithms. The Python certificate and reachable-envelope routines and the C++ causal-grid routine originate from the project's retained implementation records.

The public interface adds packaging, input handling and build instructions. The scalar release keeps reachable-envelope evaluation. The RPS release keeps the eight-round exact-prefix implementation. Computational recurrences, reward kernels, geometric feasibility tests and directed rounding retain their original definitions.

Original implementation identifiers (SHA-256):

| Component | Original source hash |
|---|---|
| Grouped-word certificate | `53c4049dc10084c1c038c23fcca8f3d0e779338a5e62321cbc8c3a455faa1a1a` |
| Scalar evaluator | `51f7c2065027cad6697fa302ffba42075f59500ea03e9c9ae8c2905b119714e7` |
| Exact-prefix RPS grid | `f60d3ead8e7773d314508470d7443db1b0bfaba0105ad1d184a23a1999bd246b` |

The RPS mathematical model follows the win-stay/lose-shift response and modal-counter strategy described in [`oudeng/RPS`](https://github.com/oudeng/RPS/tree/72c91442540027754d791eaa237fbdc60f037165). The corresponding pinned source references are [`AI_RPS/WL.py`](https://github.com/oudeng/RPS/blob/72c91442540027754d791eaa237fbdc60f037165/AI_RPS/WL.py) and [`AI_RPS/B_v1.py`](https://github.com/oudeng/RPS/blob/72c91442540027754d791eaa237fbdc60f037165/AI_RPS/B_v1.py). This repository expresses the declared model through rational state updates and finite transition/reward tables. The payoff-best-response readout is a separate adapter of that model. Upstream strategy source files are linked for attribution and are not bundled.

The rational scalar example is defined within this project. Distributed code and documentation are covered by the repository's MIT license; the linked upstream materials retain their respective ownership and terms.
