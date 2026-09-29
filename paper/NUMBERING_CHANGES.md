# Section-number concordance

The former Sections 3 and 4 were merged. The mathematical statements and Lean declaration names are unchanged. Section-based theorem and lemma numbers from former Section 5 onward have their first component reduced by one. LN97 book references, algorithm numbers, and globally numbered propositions are unaffected.

| Previous section | Current section | Previous title |
| --- | --- | --- |
| 1 | 1 | Magic Square Introduction |
| 2 | 2 | The universal form of a \(3\times 3\) magic square |
| 3 | 3 | Squares in finite fields |
| 4 | 3.2 (subsection) | Magic Squares of Center Zero and Center One |
| 5 | 4 | Quadratic characters and Weil's bound |
| 6 | 5 | A center-zero non-Parker criterion for \(q\equiv1\pmod4\) |
| 7 | 6 | A center-one non-Parker criterion for \(q\equiv3\pmod4\) |
| 8 | 7 | Pseudocode for the finite-field searches |
| 9 | 8 | Finite-field classification |
| 10 | 9 | From squares to arbitrary powers |
| 11 | 10 | Characters for \(n\)-th powers |
| 12 | 11 | The center-zero construction for \(n\)-th powers |
| 13 | 12 | Cubes |
| 14 | 13 | The center-one construction for the sign-obstructed case |
| 15 | 14 | Only finitely many odd fields fail for each fixed power |
| 16 | 15 | Characteristic two |
| 17 | 16 | Future work |
| 18 | 17 | Summary of the extension |
| 19 | 18 | Lean verification |

## Main results in the Lean guides

| Previous paper number | Current paper number | Lean declaration |
| --- | --- | --- |
| Theorem 9.1 | Theorem 8.1 | `MagicSquares.squareParker_iff_of_char_ne_two` |
| Theorem 13.1 | Theorem 12.1 | `MagicSquares.cubeParker_iff_of_char_ne_two` |
| Lemma 12.1 | Lemma 11.1 | `MagicSquares.isNthPower_neg_one_iff_powerIndex_dvd_half` |
| Theorem 12.1 | Theorem 11.1 | `MagicSquares.exists_centerZero_magic_of_powers_of_bound` |
| Theorem 14.1 | Theorem 13.1 | `MagicSquares.Square3.exists_centerOne_magic_of_powers_of_bound` |
| Theorem 15.1 | Theorem 14.1 | `MagicSquares.finite_nParker_cardinalities` |

The separate Lean repository's existing theorem map, dependency guide, and source comments may still refer to the previous paper numbering. Use this concordance until those references are updated. This paper update does not modify the Lean repository or the book's theorem numbers.
