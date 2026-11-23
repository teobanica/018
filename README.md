# Quantum partial permutations in action

**Astra 6.2 and Banica 2.0** (November 2026)

This repository contains the paper (40 pages, PDF and LaTeX source) and the 13 Python programs used for its computations.

## Summary

A partial permutation matches some points of a set and leaves the others unused. On a graph, two assignments can occur together only if they preserve adjacency and nonadjacency. Recording which assignments are compatible gives a graph of its own, the association graph.

- **Partial graph symmetries.** Cliques in the association graph describe classical partial isomorphisms. Walks give a basis for their quantum counterpart. The abstract algebra recovers the association graph, and every polynomial attains its universal norm in a finite matrix representation.
- **Rank and fixed points.** At fixed rank we determine exactly how many marked fixed-point projections can be free. The classical and quantum fixed-point laws have Poisson and free Poisson limits, respectively.
- **Magic completions and Fourier models.** We compute the smallest magic completion over two-by-two matrices and relate diagonal completions to graph coloring. We classify partial Hadamard matrices with commuting row-quotient projections above row density 5/9. Above density 2/3 they all come from Fourier matrices by deleting rows.

## Contents

[018.pdf](018.pdf) is the paper and [018.tex](018.tex) is its LaTeX source.

| Programs | Computation |
| --- | --- |
| `borders.py`, `completions.py` | Construct magic completions, verify minimum sizes and check coloring examples. |
| `ranks.py`, `fixedrank.py`, `diagonals.py` | Compute fixed-point laws and moments, construct marked projections and test diagonal identities. |
| `association.py`, `walks.py`, `reconstruct.py` | Enumerate cliques and walk words, build matrix representations and reconstruct the association graph. |
| `components.py`, `cycles.py` | Check the component relations and count partial symmetries of cycles by rank. |
| `quotients.py`, `fourcosets.py`, `translations.py` | Test commutation of row-quotient projections, enumerate four-point row sets and compute restricted Fourier translations. |

## Running the code

The programs use Python 3.9 or later and its standard library. Run these examples from the directory containing the scripts.

```sh
python3 association.py
python3 ranks.py --moment 6 4 3 4
python3 translations.py
```

The first command checks all 256 ordered pairs of two-point binary color matrices. The second computes exact classical and quantum fourth moments. The third checks the partial translations of every nonempty subset of the cyclic groups of orders 2 through 8.

