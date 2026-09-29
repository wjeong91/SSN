# CABP active-site facts for RFD3 motif design (9I27, form II Rubisco, R. rubrum)

All numbers below were measured in this project from 9I27.cif and from 95 downloaded structures
(26 PDB, 57 AlphaFold, 12 ESM). Do not invent new numbers; reason only from these.

## Reference
9I27 = "activated Form II Rubisco from Rhodospirillum rubrum with bound Magnesium and CABP".
Chains A and B, both the same large subunit. One active site is built from the C-terminal domain
of one chain plus the N-terminal domain of the partner chain. Residues below are labelled
`A<n>` (same chain as the CABP being bound) and `B<n>` (partner chain).
CABP has two phosphates (1-phosphate, 5-phosphate), a carboxylate, and O2/O3/O4 hydroxyls,
plus a Mg2+ ion.

## Every contact to CABP/Mg in 9I27 (<=3.5 A, Mg <=2.6 A)

| CABP moiety | contacts | from side chains | from Gly backbone NH |
|---|---|---|---|
| Mg | 3 | E194 OE1 2.06, D193 OD1 2.09, KCX191 OQ2 2.17 | 0 |
| 1-phosphate | 7 | T53(B) OG1 2.75 + 3.42, K329 NZ 2.88, K166 NZ 3.16 | G370 N 2.68, G394 N 2.79, G393 N 2.83 |
| 5-phosphate | 3 | R288 NE 2.85, H321 ND1 2.89, R288 NH2 3.03 | 0 (S368 backbone O sits 3.13 A away, acceptor-acceptor, likely water-mediated) |
| carboxylate | 6 | K168 NZ 2.82, N111(B) ND2 2.92, K329 NZ 2.97, D193 OD1 3.18, E194 OE1 3.24, K166 NZ 3.35 | 0 |
| O2 hydroxyl | 4 | KCX191 OQ2 3.06, K166 NZ 3.14, KCX191 OQ1 3.19, D193 OD1 3.40 | 0 |
| O3 hydroxyl | 5 | KCX191 OQ1 2.71, H287 NE2 2.83, E194 OE1 2.96, KCX191 OQ2 2.98, N111(B) ND2 3.31 | 0 |
| O4 hydroxyl | 1 | none | G369 N 3.42 (weak) |

KEY CONSEQUENCE: the glycine backbone matters for the 1-phosphate only. Every other moiety is
held entirely by side chains.

## Second-shell side-chain interactions (conserved across CABP-bound form II)

| pair | distance in 9I27 | kept in 94 CABP-bound sites | kept in AlphaFold models | what it supports |
|---|---|---|---|---|
| H285 NE2 - H321 NE2 | 2.87 A | 1.00 | 1.00 | H321, which gives 1 of 3 bonds to the 5-phosphate |
| H291 NE2 - R288 NE | 3.10 A | 1.00 | 1.00 | R288, which gives 2 of 3 bonds to the 5-phosphate. Also stacks: ring centroid to CZ 3.76 A, plane angle 17 deg |
| H44(B) NE2 - E48(B) OE | 2.86 A | 1.00 | 0.87 | E48(B), which gives ZERO direct H-bonds to CABP |
| T391 OG1 - S368 OG | 3.05 A | 0.93 | 1.00 | S368 |
| T53(B) OG1 - N54(B) OD1 | 3.40 A | 0.80 | 0.11 | T53(B) |

Side chains that are NOT within 4 A of CABP: H291 (4.06), T391 (5.72), M371 (5.96), A392 (5.10),
N54(B) (4.24), H285 (7.49), H44(B) (7.66).

## Conservation (from the user's own pipeline over 458 aligned columns)

Tier A = invariant, B = conserved, C = structurally conserved but sequence-variable, D = variable.
Tier counts across the whole protein: A 7, B 46, C 93, D 312.

| residue | tier | modal residue / fraction | note |
|---|---|---|---|
| K166 | A | K 1.00 | |
| D193 | A | D 1.00 | |
| E194 | A | E 1.00 | |
| H287 | A | H 1.00 | |
| KCX191 | B | K 0.97 | carbamylated lysine, a non-standard residue |
| R288 | B | R 0.81 | |
| H321 | B | H 0.81 | |
| K329 | B | K 0.83 | register-unstable in the MSA; Foldseek alignment says 100% K |
| S368 | B | S 0.90 | |
| G370 G393 G394 | B | G 1.00 / 0.94 / 0.82 | |
| E48(B) T53(B) N111(B) | B | E 0.81 / T 0.90 / N 0.84 | trans-chain, register-unstable, 10-25% gaps |
| K168 | C | K 0.64 family-wide, but **1.00 within form II** | |
| I164 | C | I 0.64 family-wide, but **0.94 within form II** (form I has Thr) | |
| H285 | C | H 0.74 | |
| H291 | C | H 0.60 | |
| G369 | D | G 0.77 | |
| M330 | D | M 0.36; only 0.69 even within form II | EXCLUDED from the design set |
| H44(B) | D | A 0.29 modal, His is a minority | gapfrac 0.57 — poorly determined |

Note the H44(B) row: the user's family-wide pipeline scores it Tier D with His NOT the modal
residue, while the form II structures show the H44-E48 salt bridge at 100%. These disagree; the
salt bridge is real in form II but H44 is not conserved outside it.

## Sequence positions (for grouping into contigs)

Chain A residues of interest: 164, 166, 168, 191, 193, 194, 285, 287, 288, 291, 321, 329,
368, 369, 370, 393, 394.
Chain B residues of interest: 44, 48, 53, 111.

Gaps that matter: 194->285 is 91; 291->321 is 30; 321->329 is 8; 329->368 is 39;
370->393 is 23; 48->53 is 5; 53->111 is 58.

## RFD3 input syntax (verified from RosettaCommons/foundry models/rfd3/docs/input.md)

- `contig` — indexed motif, e.g. `"A1-80,10,/0,B5-12"`. Chain-prefixed ranges come from the input;
  bare numbers are designed length; `/0` is a chain break.
- `unindex` — unindexed motif components whose sequence position is unknown to the model.
  Comma-separated components. `A11-12` ties two adjacent residues. `A11,3,A12` ties two residues
  with a stated 3-residue separation WITHOUT fixing the intervening backbone. Internal breakpoints
  are inferred and logged.
- Unindexed residues are always fixed unless `select_fixed_atoms` says otherwise; at least one atom
  of each unindexed residue must be fixed. When `select_fixed_atoms` is given as a dictionary, ONLY
  the listed atoms are carried over from the input; everything else is diffused.
- `select_fixed_atoms` — dictionary; values may be `ALL`, `TIP`, `BKBN` (= N,CA,C,O), or an explicit
  comma-separated atom list, e.g. `A15: N,CA,C,O,NE2`.
- `select_unfixed_sequence` — where sequence may change.
- `select_hbond_donor` / `select_hbond_acceptor` — atom-wise donor/acceptor conditioning flags,
  dictionary input only. Intended to make the model BUILD such an interaction, not to preserve one.
- `ligand` — ligand by PDB chemical component name (CABP is `CAP`).
- There is NO distance-restraint field. Geometry is preserved by fixing atom coordinates.
- `cleanup_guideposts: False` keeps guidepost outputs for inspection.
- Errors are raised if indexed and unindexed motifs overlap.
- CLI: `rfd3 design out_dir=<path> inputs=<inputs.yaml>`

## The constraint for this task

The user wants AT MOST 7 motif components ("조각"), because more is unrealistic for RFD3
(published examples use 2-3). Two ways to spend a component:
 (a) a contiguous span, which fixes every backbone residue in it, or
 (b) an offset-tied group like `A321,7,A329`, which fixes only the listed residues.

## Exhaustive search already run (do not redo, build on it)

Contiguous spans only, best 7-component solutions (value = sum of per-residue CABP bond weights,
64.8 = all 20 residues):
1. value 60.8, 36 residues fixed: A164-168 A191-194 A285-291 A321-329 A368-370 A393-394 B48-53 (drops B111)
2. value 59.8, 35 fixed: ... B48-53 B111 (drops G393 G394)
3. value 59.8, 31 fixed: ... A393-394 B111 (drops E48 T53)

With offset ties, 7 components can cover ALL 20 residues while fixing only 20 residues.
Smallest-maximum-span solution found:
  A164,1,A166,1,A168
  A191,1,A193-194
  A285,1,A287-288,2,A291
  A321,7,A329
  A368-370,22,A393-394
  B48,4,B53
  B111
(max span within a component = 26, at the A368...A394 group)
