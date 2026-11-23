# Quantum partial permutations in action

The 13 Python programs for paper 018 require Python 3.9 or later and use only the standard library. Unzip the archive, open the code018 folder in a terminal, and run

```sh
python3 association.py
python3 walks.py
python3 components.py
python3 fixedrank.py
python3 diagonals.py
python3 cycles.py
python3 reconstruct.py
python3 quotients.py
python3 fourcosets.py
python3 translations.py
python3 borders.py
python3 completions.py
python3 ranks.py
```

**Partial graph isomorphisms.** `association.py` builds the undirected association graph from one or two colored directed graphs. It counts its cliques by size and its walks by number of letters. With no input it checks all 256 ordered pairs of two-point binary color matrices against direct enumeration of partial bijections. It also checks the three-point single-edge graph, unequal sizes, and the directed and undirected cycle word counts through seven vertices. Clique enumeration is exponential for large inputs.

**A graph of your own.** A JSON file can contain `{"source":[[0,1],[0,0]]}` for the single arrow, or add a second square matrix as `"target"`. Vertex indices start at zero. Diagonal entries encode vertex colors, and the row is the source of a directed edge in the input adjacency matrix. Run `python3 association.py graph.json --degree 8`. The output lists allowed assignments as pairs (target, source). The classical count at rank zero includes the empty map. The quantum count at length zero includes the unit.

**Word representations.** `walks.py` constructs rational matrices for a prescribed reduced word. It checks idempotence, self-adjointness for the specified Gram inner product, and all required zero products. Its default run checks every basis word through length four on the empty two-point graph against all nonempty words up to its length. For the single-arrow example above try `python3 walks.py graph.json --walk '[[0,0],[1,1],[0,0]]'`. The displayed coefficient should be 1/16.

**Fixed-point projections at prescribed rank.** `fixedrank.py` checks all classical fixed-point sets through seven points. It constructs the universal marked projections for every admissible size and tests exact obstructions outside the theorem's range through twenty points. It also checks sample Pauli fibers for the three-point rank-two quotient and the four-point models separating reduced words through degree four, as well as the two-dimensional model attaining the three-mark norm 1/2. The module exposes `universal_array(n, rank, marks)` for constructing other sizes. For instance, `universal_array(3, 2, 2)` returns a rank-two array with two free marked projections. A positive integer denotes a free projection, its negative denotes its complement, and `"I"` denotes the identity.

**Additional full-permutation comparisons.** `diagonals.py` verifies the rational five-point magic unitary and prints nine times its nonzero degree-five diagonal identity difference. It also checks three four-point Pauli fibers and enumerates fixed-point patterns through order eight. The proof for every Pauli fiber is in the manuscript.

**Graph components.** `components.py` checks the adjacency equation against the component description for every simple graph and partial bijection on one through four vertices. It also recovers the component projections of an exact rational operator-valued example with noncommuting supports and checks all five translation fibers of a periodic version of the decorated path.

**Cycles.** `cycles.py` compares the run-profile formula for classical partial automorphisms with direct enumeration through order seven, reporting counts by rank. It checks the cyclic component formula by scanning all vertex subsets at the same orders. For larger formula-only runs try `python3 cycles.py --max-n 20 --brute-through 0`.

**Partial Hadamard matrices.** `quotients.py` verifies orthogonality and commutativity exactly by cyclotomic polynomial division and computes the finest invariant coloring of ordered row pairs. It includes the two-layer and three-layer families and a noncommuting example. To inspect another partial Butson matrix, run `python3 quotients.py matrix.json`. A sample input is `{"q":3,"rows":[[0,0,0],[0,1,2]]}`. The entry with exponent `a` means exp(2πia/q). Output indices start at zero.

**Four row cosets.** Run `python3 fourcosets.py` to enumerate four-point sets in cyclic groups through order 35 and seven product groups. It tests the progression conclusion under the seven-difference, stabilizer and torsion hypotheses. It also computes the stabilizers and the step orders of five four-layer matrices from their entries.

**Restricted translations.** `translations.py` calculates semigroups by closure under composition and builds the same sets independently by intersecting translates. By default it checks every nonempty subset of every cyclic group of orders 2 through 8. It also counts the distinct fixed-point sets for the full layered examples. Set `--max-n 9` for a larger Fourier scan.

**Magic completions.** `borders.py` constructs the rectangular completion over scalar and rational matrix coefficients. It compares attainable scalar ranks against restrictions of all full permutations through size six. It also evaluates the three tensor legs in the transpose identity separately and checks what happens in degrees one and two.

**Minimum completion size.** `completions.py` builds optimal completions with coefficients in two-by-two matrices. The default run checks every 2-by-3, 3-by-2, and 3-by-3 submagic corner drawn from zero, the identity, and three complementary pairs of lines. For each one it also tries to place every row deficit in one fewer new column. The program checks the 18-vertex example by exhaustive ordinary coloring and by verifying the four-color quantum certificate with integer inner products. To construct a corner of your own in Python, use `from completions import complete` and try `complete([[-1,0,0],[0,-2,0],[0,0,-3]])`. Opposite integers denote complementary projections and `"I"` denotes the identity.

**Rank and fixed points.** `ranks.py` compares the exact probability generating function with direct enumeration at all ranks and truncations through six vertices. To print a single law without enumeration, use `python3 ranks.py --law 20 8 6`. The output pairs a fixed-point count with its exact probability.

Use `python3 ranks.py --moment 6 4 3 4` for the exact classical and quantum fourth moments from the partition Gram matrices. The default run compares the trace formula with the enumerated classical law and recovers the full fourth moments 15 and 14.

**Reconstruction from the abstract algebra.** `reconstruct.py` starts with the character points in the projection coordinates supplied by the proof. It reconstructs the graph from equality of pair sums, tries every base of minimum degree, and checks all labeled graphs on at most five vertices by independent relabeling. It also recovers both square rank distributions and compares them with direct enumeration of partial maps.
