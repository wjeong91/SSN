# Measured backbone vs side-chain role, 17 residues of the 6-component RFD3 fallback

Measured directly from 9I27.cif (active site A: CAP A501, MG A500) with Biopython.
Cutoffs: H-bond <= 3.5 A between polar (N/O) atoms; Mg coordination <= 2.6 A;
side-chain-to-side-chain contact <= 4.0 A any heavy atom; water bridge <= 3.3 A on both legs.
SASA is Shrake-Rupley over chain A alone, in A^2.
"C needed for next" = this residue's C atom sets the backbone amide direction of residue i+1,
and residue i+1 is in the motif and donates a backbone N-H to the ligand.

The 6 components under consideration:
  A166,1,A168 | A191,1,A193-194 | A285,1,A287-288,2,A291 | A321,7,A329 | A368-370 | A392-394

MEASUREMENTS (site A of 9I27.cif)

A166 LYS
   backbone_min 7.03 (CA->O1P)   sidechain_min 3.14 (NZ->O2)  n_sidechain_atoms 5  SASA 42.0
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(3.16, 'NZ', 'O1'), (3.14, 'NZ', 'O2'), (3.35, 'NZ', 'O6')]
   sidechain_to_sidechain_within_motif: [(2.78, 'ASP193', 'NZ', 'OD2')]
   C_needed_for_next_amide: False   water_bridge: [('NZ', 2.84)]
A168 LYS
   backbone_min 6.73 (CA->O6)   sidechain_min 2.82 (NZ->O6)  n_sidechain_atoms 5  SASA 78.7
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(2.82, 'NZ', 'O6')]
   sidechain_to_sidechain_within_motif: [(3.17, 'ASP193', 'CG', 'OD2'), (2.83, 'GLU194', 'NZ', 'OE2')]
   C_needed_for_next_amide: False   water_bridge: [('NZ', 2.99)]
A191 KCX
   backbone_min 7.24 (O->O2)   sidechain_min 2.17 (OQ2->MG)  n_sidechain_atoms 8  SASA 0.0
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(3.19, 'OQ1', 'O2'), (2.71, 'OQ1', 'O3'), (3.06, 'OQ2', 'O2'), (2.98, 'OQ2', 'O3'), (2.17, 'OQ2', 'MG')]
   sidechain_to_sidechain_within_motif: [(2.81, 'ASP193', 'OQ2', 'OD1'), (2.9, 'GLU194', 'OQ2', 'OE1'), (3.35, 'HIS285', 'CE', 'CE1'), (3.41, 'HIS287', 'CX', 'CD2'), (3.35, 'HIS321', 'OQ1', 'CE1'), (3.82, 'SER368', 'OQ1', 'OG')]
   C_needed_for_next_amide: False   water_bridge: -
A193 ASP
   backbone_min 4.31 (CA->MG)   sidechain_min 2.09 (OD1->MG)  n_sidechain_atoms 4  SASA 0.0
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(3.4, 'OD1', 'O2'), (3.18, 'OD1', 'O6'), (2.09, 'OD1', 'MG')]
   sidechain_to_sidechain_within_motif: [(2.78, 'LYS166', 'OD2', 'NZ'), (3.17, 'LYS168', 'OD2', 'CG'), (2.81, 'KCX191', 'OD1', 'OQ2'), (2.88, 'GLU194', 'OD1', 'OE1')]
   C_needed_for_next_amide: False   water_bridge: -
A194 GLU
   backbone_min 3.96 (N->MG)   sidechain_min 2.06 (OE1->MG)  n_sidechain_atoms 5  SASA 30.2
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(2.96, 'OE1', 'O3'), (3.24, 'OE1', 'O6'), (2.06, 'OE1', 'MG')]
   sidechain_to_sidechain_within_motif: [(2.83, 'LYS168', 'OE2', 'NZ'), (2.9, 'KCX191', 'OE1', 'OQ2'), (2.88, 'ASP193', 'OE1', 'OD1'), (3.08, 'HIS287', 'OE1', 'NE2')]
   C_needed_for_next_amide: False   water_bridge: -
A285 HIS
   backbone_min 10.52 (C->O5P)   sidechain_min 7.59 (NE2->O3)  n_sidechain_atoms 6  SASA 1.1
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: [(3.35, 'KCX191', 'CE1', 'CE'), (2.87, 'HIS321', 'NE2', 'NE2')]
   C_needed_for_next_amide: False   water_bridge: -
A287 HIS
   backbone_min 6.03 (C->O5P)   sidechain_min 2.83 (NE2->O3)  n_sidechain_atoms 6  SASA 1.2
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(2.83, 'NE2', 'O3')]
   sidechain_to_sidechain_within_motif: [(3.41, 'KCX191', 'CD2', 'CX'), (3.08, 'GLU194', 'NE2', 'OE1'), (3.63, 'HIS321', 'CD2', 'NE2')]
   C_needed_for_next_amide: False   water_bridge: [('ND1', 3.18)]
A288 ARG
   backbone_min 4.87 (N->O5P)   sidechain_min 2.85 (NE->O5P)  n_sidechain_atoms 7  SASA 0.0
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(2.85, 'NE', 'O5P'), (3.03, 'NH2', 'O4P')]
   sidechain_to_sidechain_within_motif: [(3.1, 'HIS291', 'NE', 'NE2'), (3.97, 'HIS321', 'NE', 'CB')]
   C_needed_for_next_amide: False   water_bridge: -
A291 HIS
   backbone_min 7.20 (N->O5P)   sidechain_min 4.06 (NE2->O5P)  n_sidechain_atoms 6  SASA 3.4
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: [(3.1, 'ARG288', 'NE2', 'NE')]
   C_needed_for_next_amide: False   water_bridge: [('NE2', 2.81)]
A321 HIS
   backbone_min 4.67 (CA->O6P)   sidechain_min 2.89 (ND1->O6P)  n_sidechain_atoms 6  SASA 0.0
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(2.89, 'ND1', 'O6P')]
   sidechain_to_sidechain_within_motif: [(3.35, 'KCX191', 'CE1', 'OQ1'), (2.87, 'HIS285', 'NE2', 'NE2'), (3.63, 'HIS287', 'NE2', 'CD2'), (3.97, 'ARG288', 'CB', 'NE'), (3.68, 'SER368', 'ND1', 'CB')]
   C_needed_for_next_amide: False   water_bridge: [('ND1', 3.2)]
A329 LYS
   backbone_min 6.62 (N->O2P)   sidechain_min 2.88 (NZ->O2P)  n_sidechain_atoms 5  SASA 95.6
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: [(2.97, 'NZ', 'O7'), (2.88, 'NZ', 'O2P')]
   sidechain_to_sidechain_within_motif: -
   C_needed_for_next_amide: False   water_bridge: -
A368 SER
   backbone_min 3.13 (O->O6P)   sidechain_min 3.43 (OG->C3)  n_sidechain_atoms 2  SASA 0.0
   backbone_Hbonds_to_ligand: [(3.13, 'O', 'O6P')]
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: [(3.82, 'KCX191', 'OG', 'OQ1'), (3.68, 'HIS321', 'CB', 'ND1')]
   C_needed_for_next_amide: True   water_bridge: -
A369 GLY
   backbone_min 3.08 (CA->O4)   sidechain_min 99.00 (-->-)  n_sidechain_atoms 0  SASA 0.0
   backbone_Hbonds_to_ligand: [(3.42, 'N', 'O4')]
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: -
   C_needed_for_next_amide: True   water_bridge: -
A370 GLY
   backbone_min 2.68 (N->O2P)   sidechain_min 99.00 (-->-)  n_sidechain_atoms 0  SASA 2.4
   backbone_Hbonds_to_ligand: [(2.68, 'N', 'O2P')]
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: -
   C_needed_for_next_amide: False   water_bridge: -
A392 ALA
   backbone_min 3.84 (C->O3P)   sidechain_min 5.10 (CB->O3P)  n_sidechain_atoms 1  SASA 0.0
   backbone_Hbonds_to_ligand: -
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: -
   C_needed_for_next_amide: True   water_bridge: -
A393 GLY
   backbone_min 2.83 (N->O3P)   sidechain_min 99.00 (-->-)  n_sidechain_atoms 0  SASA 0.0
   backbone_Hbonds_to_ligand: [(2.83, 'N', 'O3P')]
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: -
   C_needed_for_next_amide: True   water_bridge: -
A394 GLY
   backbone_min 2.79 (N->O1P)   sidechain_min 99.00 (-->-)  n_sidechain_atoms 0  SASA 25.2
   backbone_Hbonds_to_ligand: [(2.79, 'N', 'O1P')]
   sidechain_Hbonds_or_Mg: -
   sidechain_to_sidechain_within_motif: -
   C_needed_for_next_amide: False   water_bridge: -
