#!/usr/bin/env python3
"""Are the side chains around the CABP site of 9I27 kept in other form II Rubisco structures?

For every structure in structures/{pdb,afdb,esm}:
  * find Rubisco large-subunit chains (pairwise alignment to 9I27 chain A)
  * build active sites: C-terminal domain of one chain + N-terminal domain of the partner chain
    (monomeric models only have the same-chain residues)
  * superpose the site on 9I27 (backbone of conserved catalytic residues; residues of the N-terminal
    domain are compared in their own local frame)
  * per residue: same amino acid?, side-chain centroid shift, chi1 difference
  * per side-chain interaction seen in 9I27: still within distance?

usage: python sidechain_conservation.py 9I27.cif [structures_dir]
"""
import gzip
import io
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.vectors import calc_dihedral

HERE = Path(__file__).resolve().parent
REF_CIF = Path(sys.argv[1])
STRUCT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "structures"

# residues whose side chain acts on the CABP site in 9I27 (form II, see form2_sidechain_influence.tsv)
RESIDUES = {
    5: ["A164", "A166", "A168", "A191", "A193", "A194", "A287", "A288", "A291", "A321", "A329", "A368",
        "B48", "B54", "B111"],
    7: ["A195", "A289", "A366", "A391", "B57"],
    9: ["A196", "A261", "A263", "A285", "A305", "B44", "B106"],
}
# side chain - side chain interactions in 9I27 (atom pairs, cutoff)
INTERACTIONS = [
    ("A285", "NE2", "A321", "NE2", 3.5, "H285-H321 H-bond"),
    ("B44", "NE2", "B48", "OE2", 4.0, "H44(B)-E48(B) salt bridge"),
    ("A391", "OG1", "A368", "OG", 3.5, "T391-S368 H-bond"),
    ("A291", "NE2", "A288", "NE", 3.6, "H291-R288 H-bond/stacking"),
    ("B54", "OD1", "B53", "OG1", 3.6, "N54(B)-T53(B) H-bond"),
    ("A195", "CD", "A194", "CB", 4.5, "P195-E194 packing"),
    ("A289", "CB", "B111", "OD1", 4.5, "A289-N111(B) packing"),
    ("A263", "OD1", "A194", "CB", 4.5, "D263-E194 packing"),
    ("A261", "CD1", "A191", "CE", 4.5, "L261-K191 packing"),
    ("B106", "CG2", "A194", "CG", 4.5, "T106(B)-E194 packing"),
]
CORE_C = [166, 191, 193, 194, 287, 288, 321, 368, 393, 394]   # superposition, C-terminal domain site
CORE_N = [48, 49, 50, 51, 53, 110, 111, 114]                    # superposition, N-terminal domain
N_DOMAIN_MAX = 140                                              # residues below this belong to the N-terminal domain
THREE = {"KCX": "K", "MSE": "M", "LLP": "K", "CSO": "C", "SEP": "S", "TPO": "T", "HYP": "P", "CME": "C"}
BB = ("N", "CA", "C", "O")

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score, aligner.extend_gap_score = -11, -1
aligner.end_insertion_score = aligner.end_deletion_score = 0


def one_letter(res):
    from Bio.Data.PDBData import protein_letters_3to1
    n = res.get_resname()
    return THREE.get(n) or protein_letters_3to1.get(n) or protein_letters_3to1.get(n.capitalize(), "X")


def load(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        text = fh.read()
    cif = ".cif" in path.name
    parser = MMCIFParser(QUIET=True) if cif else PDBParser(QUIET=True)
    return parser.get_structure(path.name.split(".")[0], io.StringIO(text))[0]


def chain_residues(chain):
    return [r for r in chain if all(a in r for a in ("N", "CA", "C"))]


def map_to_ref(res_list, ref_seq):
    seq = "".join(one_letter(r) for r in res_list)
    if not seq:
        return {}, 0, 0
    a = aligner.align(ref_seq, seq)[0]
    m = {}
    for (t0, t1), (q0, q1) in zip(*a.aligned):
        for k in range(t1 - t0):
            m[t0 + k] = q0 + k
    ident = np.mean([ref_seq[k] == seq[v] for k, v in m.items()]) if m else 0
    return m, ident, len(m) / len(ref_seq)


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((Q - qc).T @ (P - pc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return lambda X: (np.asarray(X) - qc) @ R.T + pc


# chemically equivalent / flip-ambiguous side-chain atoms
EQUIV = {"OE1": ("OE1", "OE2"), "OE2": ("OE1", "OE2"), "OD1": ("OD1", "OD2", "ND2"), "OD2": ("OD1", "OD2"),
         "ND1": ("ND1", "NE2", "CD2", "CE1"), "NE2": ("ND1", "NE2", "CD2", "CE1", "OE1"), "ND2": ("OD1", "ND2"),
         "NE": ("NE", "NH1", "NH2"), "NH1": ("NE", "NH1", "NH2"), "NH2": ("NE", "NH1", "NH2")}


def atom_dist(r1, a1, r2, a2):
    """shortest distance allowing equivalent atoms (carboxylate O, His ring N, amide flips)"""
    s1 = [r1[n].coord for n in EQUIV.get(a1, (a1,)) if n in r1]
    s2 = [r2[n].coord for n in EQUIV.get(a2, (a2,)) if n in r2]
    if not s1 or not s2:
        return None
    return float(min(np.linalg.norm(x - y) for x in s1 for y in s2))


def chi1(res):
    g = next((n for n in ("CG", "OG", "OG1", "SG", "CG1") if n in res), None)
    if g is None or "CB" not in res:
        return None
    return math.degrees(calc_dihedral(res["N"].get_vector(), res["CA"].get_vector(), res["CB"].get_vector(), res[g].get_vector()))


def sc_atoms(res):
    return {a.get_id(): a.coord for a in res if a.get_id() not in BB and a.element != "H"}


class Site:
    """one active site: C-domain chain + N-domain chain (may be the same chain for monomer models)"""

    def __init__(self, c_res, c_map, n_res, n_map, label):
        self.c_res, self.c_map, self.n_res, self.n_map, self.label = c_res, c_map, n_res, n_map, label

    def residue(self, tag, offset=3):
        chain, num = tag[0], int(tag[1:])
        res, m = (self.c_res, self.c_map) if chain == "A" else (self.n_res, self.n_map)
        k = num - offset
        return res[m[k]] if k in m else None


def reference_sites():
    model = MMCIFParser(QUIET=True).get_structure("ref", str(REF_CIF))[0]
    A, B = chain_residues(model["A"]), chain_residues(model["B"])
    seq = "".join(one_letter(r) for r in A)
    idm = {i: i for i in range(len(A))}
    return seq, Site(A, idm, B, idm, "9I27 A/B")


def frame(site_ref, site, core, chain):
    P, Q = [], []
    for num in core:
        r0, r1 = site_ref.residue(f"{chain}{num}"), site.residue(f"{chain}{num}")
        if r0 is None or r1 is None:
            continue
        for a in ("N", "CA", "C"):
            P.append(r0[a].coord)
            Q.append(r1[a].coord)
    return kabsch(np.array(P), np.array(Q)) if len(P) >= 15 else None


def compare(site_ref, site):
    f_c = frame(site_ref, site, CORE_C, "A")
    f_n = frame(site_ref, site, CORE_N, "B") if site.n_res is not site.c_res else frame(site_ref, site, CORE_N, "B")
    out = {}
    for tags in RESIDUES.values():
        for tag in tags:
            r0, r1 = site_ref.residue(tag), site.residue(tag)
            f = f_c if tag[0] == "A" else f_n
            if r1 is None or f is None:
                out[tag] = None
                continue
            same = one_letter(r0) == one_letter(r1)
            s0, s1 = sc_atoms(r0), sc_atoms(r1)
            common = [a for a in s0 if a in s1]
            shift = float(np.linalg.norm(np.mean([s0[a] for a in common], 0) - np.mean(f([s1[a] for a in common]), 0))) if common else None
            c0, c1 = chi1(r0), chi1(r1)
            dchi = None if c0 is None or c1 is None else abs((c1 - c0 + 180) % 360 - 180)
            out[tag] = dict(aa=one_letter(r1), same=same, shift=shift, dchi=dchi)
    inter = {}
    for t1, a1, t2, a2, cut, name in INTERACTIONS:
        r1, r2 = site.residue(t1), site.residue(t2)
        cross = t1[0] != t2[0] and site.n_res is site.c_res   # needs two chains
        d = None if r1 is None or r2 is None or cross else atom_dist(r1, a1, r2, a2)
        inter[name] = None if d is None else d <= cut
    return out, inter


def state(model, site):
    """ligand / carbamylation / Mg near this site"""
    k191 = site.residue("A191")
    if k191 is None:
        return "n/a"
    center = k191["CA"].coord
    ligs = set()
    for r in model.get_residues():
        if r.id[0].startswith("H_") and r.get_resname() not in THREE:
            if min(np.linalg.norm(a.coord - center) for a in r) < 10:
                ligs.add(r.get_resname())
    carb = k191.get_resname() == "KCX"
    return f"{'+'.join(sorted(ligs)) or 'apo'}{' KCX' if carb else ''}"


def sites_of(model, ref_seq):
    chains = []
    for ch in model:
        res = chain_residues(ch)
        m, ident, cov = map_to_ref(res, ref_seq)
        if ident >= 0.45 and cov >= 0.6:
            chains.append((ch.id, res, m))
    sites = []
    for cid, res, m in chains:
        k = 191 - 3
        if k not in m:
            continue
        anchor = res[m[k]]["CA"].coord
        best = None
        for cid2, res2, m2 in chains:
            if cid2 == cid:
                continue
            pts = [res2[m2[n - 3]]["CA"].coord for n in (48, 53, 111) if (n - 3) in m2]
            if len(pts) < 2:
                continue
            d = float(np.mean([np.linalg.norm(p - anchor) for p in pts]))
            if best is None or d < best[0]:
                best = (d, res2, m2, cid2)
        if best and best[0] < 20:
            sites.append(Site(res, m, best[1], best[2], f"{cid}/{best[3]}"))
        else:
            sites.append(Site(res, m, res, m, f"{cid} (monomer)"))
    return sites


def main():
    ref_seq, ref_site = reference_sites()
    files = sorted(p for p in STRUCT.rglob("*") if p.suffix in (".gz", ".cif", ".pdb") and p.is_file() and "manifest" not in p.name)
    per_res = defaultdict(lambda: defaultdict(list))
    per_int = defaultdict(lambda: defaultdict(list))
    rows = []
    for path in files:
        source = path.parent.name
        try:
            model = load(path)
        except Exception as e:  # noqa: BLE001
            print(f"skip {path.name}: {e}")
            continue
        for site in sites_of(model, ref_seq):
            name = path.name.split(".")[0].upper()
            if name == "9I27" and site.label.startswith("A"):
                continue  # the reference site itself
            res_out, int_out = compare(ref_site, site)
            st = state(model, site) if source == "pdb" else "model"
            group = source if source != "pdb" else ("pdb_CABP" if "CAP" in st.split(" ")[0].split("+") else "pdb_noCABP")
            if name == "9I27":
                group = "control_9I27_siteB"
            rows.append((path.name.split(".")[0], site.label, source, st))
            for tag, v in res_out.items():
                if v is not None:
                    v["entry"] = name
                per_res[tag][group].append(v)
            for name, v in int_out.items():
                per_int[name][group].append(v)
    groups = ["control_9I27_siteB", "pdb_CABP", "pdb_noCABP", "afdb", "esm"]

    with open(HERE / "sidechain_sites.tsv", "w") as fh:
        fh.write("structure\tsite\tsource\tstate\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    with open(HERE / "sidechain_residues.tsv", "w") as fh:
        fh.write("shell_A\tresidue\tgroup\tn\tsame_aa\tsc_shift_median\tsc_within_1.5A\tchi1_same_rotamer\tkept\n")
        print(f"{'res':8s} " + " ".join(f"{g:>22s}" for g in groups))
        for shell, tags in RESIDUES.items():
            for tag in tags:
                cells = []
                for g in groups:
                    vals = [v for v in per_res[tag][g] if v is not None]
                    if not vals:
                        cells.append(f"{'-':>22s}")
                        continue
                    same = np.mean([v["same"] for v in vals])
                    sh = [v["shift"] for v in vals if v["same"] and v["shift"] is not None]
                    ch = [v["dchi"] for v in vals if v["same"] and v["dchi"] is not None]
                    within = np.mean([s <= 1.5 for s in sh]) if sh else float("nan")
                    rot = np.mean([c <= 60 for c in ch]) if ch else float("nan")
                    kept = np.mean([v["same"] and v["shift"] is not None and v["shift"] <= 1.5 and (v["dchi"] is None or v["dchi"] <= 60) for v in vals])
                    n_entries = len({v["entry"] for v in vals})
                    fh.write(f"{shell}\t{tag}\t{g}\t{len(vals)} sites/{n_entries} entries\t{same:.2f}\t{np.median(sh) if sh else float('nan'):.2f}\t{within:.2f}\t{rot:.2f}\t{kept:.2f}\n")
                    cells.append(f"n{len(vals):3d} kept {kept:4.2f} ({np.median(sh) if sh else float('nan'):3.1f}A)")
                label = f"{ref_seq[int(tag[1:]) - 3]}{tag[1:]}" + ("(B)" if tag[0] == "B" else "")
                print(f"{label:8s} " + " ".join(f"{c:>22s}" for c in cells))

    with open(HERE / "sidechain_interactions.tsv", "w") as fh:
        fh.write("interaction\tgroup\tn_sites\tfraction_present\n")
        print()
        for name in [i[5] for i in INTERACTIONS]:
            cells = []
            for g in groups:
                vals = [v for v in per_int[name][g] if v is not None]
                if vals:
                    fh.write(f"{name}\t{g}\t{len(vals)}\t{np.mean(vals):.2f}\n")
                    cells.append(f"n{len(vals):3d} {np.mean(vals):4.2f}")
                else:
                    cells.append("-")
            print(f"{name:28s} " + " ".join(f"{c:>12s}" for c in cells))
    print(f"\n{len(rows)} active sites from {len(files)} files; groups: {', '.join(groups)}")


if __name__ == "__main__":
    main()
