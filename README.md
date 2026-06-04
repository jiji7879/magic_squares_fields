# Magic Squares of Squares over Finite Fields

This repository contains the Python code and output files used for the computational verification in the paper:

**Magic Squares of Squares over Finite Fields and an Extension to Arbitrary Powers**

The main problem is to determine which finite fields are **Parker fields**, meaning finite fields over which there is no ($3\times 3$) magic square whose nine entries are distinct squares.

The computations here support the finite-field classification by checking the remaining fields below the theoretical character-sum bounds.

## Summary of results

The complete list of Parker finite fields for squares is:

\[
\mathbb F_3,\mathbb F_5,\mathbb F_7,\mathbb F_9,\mathbb F_{11},\mathbb F_{13},\mathbb F_{17},
\mathbb F_{19},\mathbb F_{23},\mathbb F_{25},\mathbb F_{27},\mathbb F_{31},
\mathbb F_{43},\mathbb F_{47},\mathbb F_{67},\mathbb F_{243}.
\]

For fields of order (Q\equiv 1 \pmod 4), the Parker fields are:

\[
\mathbb F_5,\mathbb F_9,\mathbb F_{13},\mathbb F_{17},\mathbb F_{25}.
\]

For fields of order (Q\equiv 3 \pmod 4), the Parker fields are:

\[
\mathbb F_3,\mathbb F_7,\mathbb F_{11},\mathbb F_{19},\mathbb F_{23},
\mathbb F_{27},\mathbb F_{31},\mathbb F_{43},\mathbb F_{47},
\mathbb F_{67},\mathbb F_{243}.
\]

The repository also includes a computation for cubes. Among odd-characteristic finite fields, the (3)-Parker fields found below the theoretical bound (Q<1037) are:

\[
\mathbb F_3,\mathbb F_5,\mathbb F_7,\mathbb F_{13},\mathbb F_{19},\mathbb F_{25},
\mathbb F_{31},\mathbb F_{37},\mathbb F_{43},\mathbb F_{61},\mathbb F_{67},
\mathbb F_{79},\mathbb F_{127},\mathbb F_{343}.
\]

Characteristic two is exceptional: every ($3\times3$) magic square over a field of characteristic two has repeated entries, so the distinct-entry problem has permanent characteristic-two obstructions.

## Files

| File                           | Purpose                                                                           |
| ------------------------------ | --------------------------------------------------------------------------------- |
| `pythagorean_fields_common.py`    | Shared finite-field arithmetic and magic-square search functions.                                               |
| `pythagorean_fields_q1.py`        | Verifies the ($Q\equiv 1\mod 4$) square case below the theoretical bound (Q<77).                                          |
| `q1_results.txt`                  | Output from the ($Q\equiv 1\mod 4$) square computation.                                                                 |
| `pythagorean_fields_q3.py`        | Verifies the ($Q\equiv 3\mod 4$) square case below the theoretical bound (Q<553736).                                      |
| `q3_b3_then_autoline_results.txt` | Output from the ($Q\equiv 3\mod4$) square computation finding witnesses with the b3 then autoline construction.        |
| `q3_full_uv_results.txt`          | Output from the ($Q\equiv 3\mod4$) square computation for the fields without witnesses in q3_b3_then_autoline_results. |
| `pythagorean_fields_n3.py`        | Verifies the odd-characteristic cube case below the theoretical bound ($Q<1037$).                                        |
| `n3_results.txt`                  | Output from the cube computation.                                                                 |

## Mathematical background

Every ($3\times3$) magic square over a commutative ring can be written in the universal form

\[
M(E,U,V)=
\begin{pmatrix}
E+U & E-U-V & E+V\
E-U+V & E & E+U-V\
E-V & E+U+V & E-U
\end{pmatrix}.
\]

For finite fields of odd order, the search can be normalized into center-zero and center-one cases.

When ($Q\equiv 1\mod 4$), (-1) is a square in ($\mathbb F_Q$), so the center-zero construction is especially useful.

When ($Q\equiv 3\mod 4$), (-1) is not a square in ($\mathbb F_Q$), so the center-one construction is used.

The paper proves that:

* every ($Q\equiv 1\mod 4$) finite field with ($Q\ge 77$) is non-Parker;
* every ($Q\equiv 3\mod 4$) finite field with ($Q\ge 553736$) is non-Parker;
* every odd finite field with ($Q\ge 1037$) is non-(3)-Parker for cubes.

Thus only finitely many finite fields need to be checked computationally in each case.

## Running the computations

The scripts are written in pure Python and do not require SageMath.

Run the ($Q\equiv 1\mod 4$) square computation:

```bash
python pythagorean_fields_q1.py
```

Run the ($Q\equiv 3\mod 4$) square computation:

```bash
python pythagorean_fields_q3.py
```

Run the cube computation:

```bash
python pythagorean_fields_n3.py
```

Each script writes a corresponding `.txt` output file.

## Notes on verification

A witness proves that a field is non-Parker. For fields where no witness is found, the relevant script performs a complete normalized search over the appropriate magic-square family.

The ($Q\equiv 3\mod 4$) computation checks prime powers ($Q<553736$) satisfying

\[
Q\equiv 3\mod 4,
\]

not merely prime fields ($\mathbb F_p$) with ($p\equiv 3 \mod 4$).

This distinction matters because finite fields include extension fields such as

\[
\mathbb F_{27},\mathbb F_{243},\mathbb F_{343},\ldots.
\]

## Citation

If using this repository, please cite the accompanying paper.
