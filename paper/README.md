# Paper

**Magic Squares of Squares over Finite Fields and an Extension to Arbitrary Powers**  
David Lai

- [Read the PDF](magic-squares-over-finite-fields.pdf)
- [LaTeX source](main.tex)
- [Bibliography](finite_square_squares.bib)

The original Sections 3 and 4 have been merged into Section 3, "Squares over finite fields and normalized constructions". The Lean verification section is now Section 18. Later section-based theorem and lemma numbers shift accordingly; see [the numbering guide](NUMBERING_CHANGES.md). The mathematics is unchanged. The PDF has been rebuilt, and the merged opening sections have been visually checked.

## Build the paper

Run from this directory with a LaTeX distribution providing pdfLaTeX, BibTeX, and the packages used in `main.tex` (including TikZ, `algorithm`, and `algpseudocode`):

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Check each command succeeds. The result is `main.pdf`. Review it before replacing the distributed `magic-squares-over-finite-fields.pdf`. In PowerShell, after reviewing:

```powershell
Copy-Item main.pdf magic-squares-over-finite-fields.pdf
```

The source references `finite_square_squares.bib`; no external figures or additional source inputs were found. Since the source has no explicit `\date`, a rebuild uses the build date on the title page and need not reproduce the supplied PDF byte for byte.

## Lean formalization

The formal proofs and their build configuration live in the separate [Lean repository](https://github.com/jiji7879/magic-squares-fields-lean). Start with its [theorem map](https://github.com/jiji7879/magic-squares-fields-lean/blob/main/THEOREM_MAP.md) and [proof dependency guide](https://github.com/jiji7879/magic-squares-fields-lean/blob/main/PROOF_DEPENDENCIES.md).

On 2026-09-28, [GitHub Actions run 36498472790](https://github.com/jiji7879/magic-squares-fields-lean/actions/runs/36498472790) completed successfully for [commit `93efc5597843`](https://github.com/jiji7879/magic-squares-fields-lean/commit/93efc5597843c9da7e5ce26176f90a40ea7508d6). The project/compatibility build and endpoint-axiom check steps both passed. This is evidence for that Lean commit; it does not establish a line-by-line correspondence between this PDF and every Lean declaration. Consult the theorem map for formalization scope.

## Licensing

The MIT software license has been selected for the code. David has not yet selected a separate paper license; this packaging step does not select one or change the repository's existing license file.
