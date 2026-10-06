#!/usr/bin/env python3
"""Prepare Rubisco motif PDBs (bind_core / bind_c40 / bind_c45) for AP Novo.

stdlib only (python >= 3.8), so it runs on a cluster login node without installs.

Per input it writes <out>/<stem>/:
  motif.pdb        protein chains A/B kept as-is, KCX191 -> LYS + FMT ligand,
                   CAP -> chain C, MG -> chain D, FMT -> chain E (1 residue each)
  motif_spec.json  residue list, token counts, AP Novo manifest *fragments*
  report.txt       checks (same text as stdout)

It does NOT write a full AP Novo manifest: only the fields quoted from an earlier
reading of the README (motif_str / seq_length / unindexed_motif_residues /
is_author_naming) are emitted, as fragments. Check them against the README.

Usage:
  python3 prep_ap_novo.py --inputs inputs/*.pdb --out out [--length 360]
                          [--numbering-map map.json]
  python3 prep_ap_novo.py --selftest
"""
import argparse
import json
import math
import os
import sys

AA_HEAVY = {
    'ALA': 'N CA C O CB', 'ARG': 'N CA C O CB CG CD NE CZ NH1 NH2',
    'ASN': 'N CA C O CB CG OD1 ND2', 'ASP': 'N CA C O CB CG OD1 OD2',
    'CYS': 'N CA C O CB SG', 'GLN': 'N CA C O CB CG CD OE1 NE2',
    'GLU': 'N CA C O CB CG CD OE1 OE2', 'GLY': 'N CA C O',
    'HIS': 'N CA C O CB CG ND1 CD2 CE1 NE2',
    'ILE': 'N CA C O CB CG1 CG2 CD1', 'LEU': 'N CA C O CB CG CD1 CD2',
    'LYS': 'N CA C O CB CG CD CE NZ', 'MET': 'N CA C O CB CG SD CE',
    'PHE': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ', 'PRO': 'N CA C O CB CG CD',
    'SER': 'N CA C O CB OG', 'THR': 'N CA C O CB OG1 CG2',
    'TRP': 'N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL': 'N CA C O CB CG1 CG2',
}
AA_HEAVY = {k: v.split() for k, v in AA_HEAVY.items()}

KCX_TO_FMT = {'CX': 'C', 'OQ1': 'O1', 'OQ2': 'O2'}   # carbamate -> formate (CCD FMT)
LIGAND_CHAIN = {'CAP': 'C', 'MG': 'D', 'FMT': 'E'}
MG_CONTACT_CUTOFF = 2.6


def parse_pdb(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line[:6] not in ('ATOM  ', 'HETATM'):
                continue
            if line[16].strip():
                raise ValueError('altloc present, not handled: ' + line.rstrip())
            atoms.append({
                'rec': line[:6].strip(), 'name': line[12:16].strip(),
                'res': line[17:20].strip(), 'chain': line[21],
                'num': int(line[22:26]), 'icode': line[26].strip(),
                'xyz': (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                'occ': float(line[54:60]), 'b': float(line[60:66]),
                'elem': line[76:78].strip(),
            })
    return atoms


def group_residues(atoms):
    """-> ordered list of ((chain, num, icode, res), [atoms])."""
    out, index = [], {}
    for a in atoms:
        key = (a['chain'], a['num'], a['icode'], a['res'])
        if key not in index:
            index[key] = len(out)
            out.append((key, []))
        out[index[key]][1].append(a)
    return out


def dist(p, q):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(p, q)))


def fmt_atom(serial, rec, atom, res, chain, num):
    name, elem = atom['name'], atom['elem']
    nm = name.ljust(4) if (len(elem) == 2 or len(name) == 4) else (' ' + name).ljust(4)
    x, y, z = atom['xyz']
    return '%-6s%5d %s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s' % (
        rec, serial, nm, res, chain, num, x, y, z, atom['occ'], atom['b'], elem)


def convert(atoms):
    """Split into prepared protein / ligand residues. Returns (protein, ligands, notes)."""
    protein, ligands, notes = [], [], []
    fmt_atoms = []
    for (chain, num, icode, res), ras in group_residues(atoms):
        rid = '%s%d' % (chain, num)
        if res == 'KCX':
            by = {a['name']: a for a in ras}
            lys = [by[n] for n in AA_HEAVY['LYS'] if n in by]
            missing = [n for n in AA_HEAVY['LYS'] if n not in by]
            cx = {KCX_TO_FMT[n]: dict(by[n], name=KCX_TO_FMT[n]) for n in KCX_TO_FMT if n in by}
            if missing or len(cx) != 3:
                raise ValueError('KCX %s incomplete: lys missing %s, carbamate %s' % (rid, missing, sorted(cx)))
            protein.append({'id': rid, 'chain': chain, 'num': num, 'res': 'LYS',
                            'atoms': lys, 'orig_res': 'KCX'})
            fmt_atoms.append({'src': rid, 'atoms': [cx['C'], cx['O1'], cx['O2']]})
            notes.append('%s KCX -> LYS + FMT(chain %s); NZ-CX bond %.2f A is not encoded' % (
                rid, LIGAND_CHAIN['FMT'], dist(by['NZ']['xyz'], by['CX']['xyz'])))
        elif res in LIGAND_CHAIN and res != 'FMT':
            ligands.append({'res': res, 'chain': LIGAND_CHAIN[res], 'src': rid, 'atoms': ras})
        elif res in AA_HEAVY:
            protein.append({'id': rid, 'chain': chain, 'num': num, 'res': res,
                            'atoms': ras, 'orig_res': res})
        else:
            raise ValueError('unhandled residue %s %s' % (res, rid))
    if len(fmt_atoms) > 1:
        raise ValueError('more than one KCX; FMT chain assignment needs extending')
    for f in fmt_atoms:
        ligands.append({'res': 'FMT', 'chain': LIGAND_CHAIN['FMT'], 'src': f['src'], 'atoms': f['atoms']})
    ligands.sort(key=lambda l: l['chain'])
    return protein, ligands, notes


def write_pdb(path, protein, ligands):
    lines, serial = [], 0
    prev = None
    for r in protein:
        if prev is not None and r['chain'] != prev:
            lines.append('TER')
        prev = r['chain']
        for a in r['atoms']:
            serial += 1
            lines.append(fmt_atom(serial, 'ATOM', a, r['res'], r['chain'], r['num']))
    lines.append('TER')
    for l in ligands:
        for a in l['atoms']:
            serial += 1
            lines.append(fmt_atom(serial, 'HETATM', a, l['res'], l['chain'], 1))
        lines.append('TER')
    lines.append('END')
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def layout_motif_str(residues, design_pos, length, ligand_tail):
    """residues: ordered [(chain, num)], design_pos: {'A164': 25, ...}.
    Spacer = gap in design numbering; consecutive native AND design -> 'A368-370'."""
    items = []
    for chain, num in residues:
        items.append((design_pos['%s%d' % (chain, num)], chain, num))
    items.sort()
    pos = [p for p, _, _ in items]
    if len(set(pos)) != len(pos):
        raise ValueError('duplicate design positions')
    if pos[0] < 1 or pos[-1] > length:
        raise ValueError('design positions outside 1..%d' % length)
    parts, prev_pos, i = [], 0, 0
    while i < len(items):
        p, chain, num = items[i]
        if p - prev_pos - 1 > 0:
            parts.append(str(p - prev_pos - 1))
        j = i
        while (j + 1 < len(items) and items[j + 1][1] == chain
               and items[j + 1][2] == items[j][2] + 1 and items[j + 1][0] == items[j][0] + 1):
            j += 1
        parts.append('%s%d' % (chain, num) if j == i else '%s%d-%d' % (chain, num, items[j][2]))
        prev_pos = items[j][0]
        i = j + 1
    if length - prev_pos > 0:
        parts.append(str(length - prev_pos))
    return ','.join(parts) + ligand_tail


def selftest():
    # Layout quoted in the earlier analysis (T13 on A10, 360 res). A10 numbers are
    # second-hand from that note; this only tests the gap arithmetic.
    a10 = {'A164': 25, 'A166': 27, 'A168': 29, 'A288': 125, 'A321': 164, 'A329': 172,
           'A368': 211, 'A369': 212, 'A370': 213, 'A393': 252, 'A394': 253,
           'B53': 293, 'B111': 336}
    res = [(k[0], int(k[1:])) for k in a10]
    got = layout_motif_str(res, a10, 360, '/C1/D1')
    want = '24,A164,1,A166,1,A168,95,A288,38,A321,7,A329,38,A368-370,38,A393-394,39,B53,42,B111,24/C1/D1'
    assert got == want, '\n got  %s\n want %s' % (got, want)
    print('selftest ok:', got)


def check_and_report(name, atoms, protein, ligands, notes, length):
    out = []
    p = out.append
    p('== %s' % name)
    chains = {}
    for r in protein:
        chains.setdefault(r['chain'], []).append(r)
    for c, rs in sorted(chains.items()):
        p('protein chain %s: %d residues: %s' % (c, len(rs), ' '.join('%s%d' % (r['res'], r['num']) for r in rs)))
    for l in ligands:
        p('ligand chain %s: %s (from %s), %d atoms' % (l['chain'], l['res'], l['src'], len(l['atoms'])))
    n_lig_atoms = sum(len(l['atoms']) for l in ligands)
    n_tokens = len(protein) + n_lig_atoms
    p('motif tokens = %d residues + %d ligand atoms = %d (README cap for unindexed: 128, per earlier note)' % (
        len(protein), n_lig_atoms, n_tokens))
    p('design tokens ~ length %d + %d ligand atoms = %d (training crop 512, per earlier note)' % (
        length, n_lig_atoms, length + n_lig_atoms))
    bad = 0
    for r in protein:
        have = {a['name'] for a in r['atoms']}
        miss = [n for n in AA_HEAVY[r['res']] if n not in have]
        extra = sorted(have - set(AA_HEAVY[r['res']]))
        if miss or extra:
            bad += 1
            p('  ATOM-SET MISMATCH %s %s: missing %s extra %s' % (r['id'], r['res'], miss, extra))
    p('residue atom sets: %s' % ('all complete' if not bad else '%d residue(s) off' % bad))
    cap = [l for l in ligands if l['res'] == 'CAP']
    if cap:
        names = sorted(a['name'] for a in cap[0]['atoms'])
        p('CAP atom names (%d): %s' % (len(names), ' '.join(names)))
        p('  NOT verified against the CCD (RCSB unreachable from this container)')
    mg = [a for a in atoms if a['res'] == 'MG']
    if mg:
        contacts = sorted(
            (dist(mg[0]['xyz'], a['xyz']), '%s%d %s %s' % (a['chain'], a['num'], a['res'], a['name']))
            for a in atoms if a is not mg[0] and a['elem'] in ('O', 'N')
            and dist(mg[0]['xyz'], a['xyz']) <= MG_CONTACT_CUTOFF)
        p('Mg O/N contacts <= %.1f A: %d' % (MG_CONTACT_CUTOFF, len(contacts)))
        for d, s in contacts:
            p('  %.2f  %s' % (d, s))
    for n in notes:
        p('note: ' + n)
    return out


def residue_key(a):
    return (a['chain'], a['num'], a['res'], a['name'])


def compare(structs):
    out = ['== cross-structure check']
    keysets = {n: {residue_key(a): a['xyz'] for a in atoms} for n, atoms in structs.items()}
    names = list(keysets)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            shared = set(keysets[a]) & set(keysets[b])
            dev = max((dist(keysets[a][k], keysets[b][k]) for k in shared), default=float('nan'))
            ra = {(k[0], k[1]) for k in keysets[a]}
            rb = {(k[0], k[1]) for k in keysets[b]}
            rel = ('%s subset of %s' % (a, b)) if keysets[a].keys() <= keysets[b].keys() else (
                ('%s subset of %s' % (b, a)) if keysets[b].keys() <= keysets[a].keys() else 'not nested')
            out.append('%s vs %s: %d shared atoms, max deviation %.3f A, %s; only-in-%s: %s; only-in-%s: %s' % (
                a, b, len(shared), dev, rel, a, sorted(ra - rb), b, sorted(rb - ra)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs', nargs='+')
    ap.add_argument('--out', default='out')
    ap.add_argument('--length', type=int, default=360, help='designed protein length (default 360, A10 size)')
    ap.add_argument('--numbering-map', help='JSON {"A164": 25, ...}: native residue -> design position; '
                                            'enables the indexed motif_str')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return 0
    if not args.inputs:
        ap.error('--inputs required')
    dmap = json.load(open(args.numbering_map)) if args.numbering_map else None
    structs, report = {}, []
    for path in args.inputs:
        stem = os.path.splitext(os.path.basename(path))[0]
        atoms = parse_pdb(path)
        structs[stem] = atoms
        protein, ligands, notes = convert(atoms)
        d = os.path.join(args.out, stem)
        os.makedirs(d, exist_ok=True)
        write_pdb(os.path.join(d, 'motif.pdb'), protein, ligands)
        rep = check_and_report(stem, atoms, protein, ligands, notes, args.length)
        ids = [r['id'] for r in protein]
        tail = ''.join('/%s1' % l['chain'] for l in ligands)
        spec = {
            'source': os.path.basename(path),
            'protein_residues': [{'id': r['id'], 'resname': r['res'], 'orig_resname': r['orig_res'],
                                  'n_atoms': len(r['atoms'])} for r in protein],
            'ligands': [{'chain': l['chain'], 'resname': l['res'], 'from': l['src'],
                         'n_atoms': len(l['atoms'])} for l in ligands],
            'n_motif_tokens': len(protein) + sum(len(l['atoms']) for l in ligands),
            'unverified_fragments_from_earlier_readme_reading': {
                'unindexed': {'motif_str': '%d%s' % (args.length, tail), 'seq_length': str(args.length),
                              'unindexed_motif_residues': ids},
                'indexed': None,
            },
        }
        if dmap is not None:
            miss = [i for i in ids if i not in dmap]
            if miss:
                rep.append('indexed: skipped, numbering map lacks %s' % miss)
            else:
                ms = layout_motif_str([(r['chain'], r['num']) for r in protein], dmap, args.length, tail)
                spec['unverified_fragments_from_earlier_readme_reading']['indexed'] = {
                    'motif_str': ms, 'seq_length': str(args.length), 'is_author_naming': True}
                rep.append('indexed motif_str: ' + ms)
        else:
            rep.append('indexed: skipped (no --numbering-map)')
        json.dump(spec, open(os.path.join(d, 'motif_spec.json'), 'w'), indent=2)
        open(os.path.join(d, 'report.txt'), 'w').write('\n'.join(rep) + '\n')
        report += rep + ['']
    if len(structs) > 1:
        cmp_lines = compare(structs)
        open(os.path.join(args.out, 'cross_check.txt'), 'w').write('\n'.join(cmp_lines) + '\n')
        report += cmp_lines
    print('\n'.join(report))
    return 0


if __name__ == '__main__':
    sys.exit(main())
