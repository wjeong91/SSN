# RFD3 motif for the 9I27 CABP site, split by how each residue reaches the ligand.
#
#   object1  side-chain contributors  K166 K168 KCX191 D193 E194 H285 H287 R288 H291 H321 K329
#   object2  backbone contributors    S368 G369 G370 A392 G393 G394
#   object3  trans-chain, chain B     E48 T53 N111   (7-component spec only, starts hidden)
#
# The atoms RFD3 pins are drawn thick and coloured. The rest of each residue stays
# thin and grey. Site is CAP A501 / MG A500. Every distance drawn is one measured in
# RFD3_ATOM_MEASUREMENTS.md.
#
# HOW TO RUN
#   Load the structure FIRST, under any name, then run this file:
#       pymol 9I27_Mg_CABP.pdb
#       PyMOL> @show_rfd3_motif.pml
#   or in one line from a shell:
#       pymol your_file.pdb show_rfd3_motif.pml
#
#   This script loads nothing itself and never refers to your object by name, so the
#   file can be called anything. If nothing appears it prints why at the end.

hide everything
remove solvent

# ---------------------------------------------------------------- objects
create object1, chain A and resi 166+168+191+193+194+285+287+288+291+321+329
create object2, chain A and resi 368+369+370+392+393+394
create object3, chain B and resi 48+53+111
create ligand, chain A and resn CAP+MG
create context, polymer

# ---------------------------------------------------------------- pinned atoms
# Built one line at a time - no backslash continuations, which some PyMOL builds mishandle.
select fix1, object1 and resi 166+168+329 and name CE+NZ
select fix1, fix1 or (object1 and resi 191 and name CX+OQ1+OQ2)
select fix1, fix1 or (object1 and resi 193 and name CG+OD1+OD2)
select fix1, fix1 or (object1 and resi 194 and name CD+OE1+OE2)
select fix1, fix1 or (object1 and resi 285+287+291+321 and name ND1+CE1+NE2)
select fix1, fix1 or (object1 and resi 288 and name NE+CZ+NH2)

select fix2, object2 and resi 368 and name C+O
select fix2, fix2 or (object2 and resi 369+393 and name N+CA+C)
select fix2, fix2 or (object2 and resi 370+394 and name N+CA)
select fix2, fix2 or (object2 and resi 392 and name C)

# ---------------------------------------------------------------- style
show sticks, object1
show sticks, object2
show sticks, object3
show sticks, ligand and not resn MG
show spheres, ligand and resn MG
set sphere_scale, 0.40, ligand and resn MG
show cartoon, context
color grey90, context
set cartoon_transparency, 0.8, context

set stick_radius, 0.10, object1
set stick_radius, 0.10, object2
set stick_radius, 0.10, object3
set stick_radius, 0.25, fix1
set stick_radius, 0.25, fix2
set stick_radius, 0.18, ligand

color grey60, object1
color grey60, object2
color grey60, object3
util.cnc("ligand")
color yellow, ligand and elem C
# Mg is a dark sphere, not magenta - the mgbridge pair lines below are purple and the
# two would be hard to tell apart.
color grey20, ligand and resn MG
label ligand and resn MG, "Mg"
set label_color, grey20, ligand and resn MG

# Category colour covers the whole pinned atom set, heteroatoms included, so the
# side-chain / backbone split reads at a glance. For CPK heteroatoms instead, run:
#     util.cnc("fix1 or fix2")
color marine, fix1
color salmon, fix2
color palegreen, object3

# ---------------------------------------------------------------- measured contacts
# Side chain to CABP or Mg.
dist d_sc, object1 and resi 166 and name NZ, ligand and name O2
dist d_sc, object1 and resi 168 and name NZ, ligand and name O6
dist d_sc, object1 and resi 191 and name OQ2, ligand and resn MG
dist d_sc, object1 and resi 191 and name OQ1, ligand and name O3
dist d_sc, object1 and resi 193 and name OD1, ligand and resn MG
dist d_sc, object1 and resi 194 and name OE1, ligand and resn MG
dist d_sc, object1 and resi 287 and name NE2, ligand and name O3
dist d_sc, object1 and resi 288 and name NE, ligand and name O5P
dist d_sc, object1 and resi 288 and name NH2, ligand and name O4P
dist d_sc, object1 and resi 321 and name ND1, ligand and name O6P
dist d_sc, object1 and resi 329 and name NZ, ligand and name O2P

# Backbone to the phosphates and to O4.
dist d_bb, object2 and resi 368 and name O, ligand and name O6P
dist d_bb, object2 and resi 369 and name N, ligand and name O4
dist d_bb, object2 and resi 370 and name N, ligand and name O2P
dist d_bb, object2 and resi 393 and name N, ligand and name O3P
dist d_bb, object2 and resi 394 and name N, ligand and name O1P


# ---------------------------------------------------------------- side-chain pairs
# Every polar side-chain pair inside the motif, measured in 9I27 and classified by the
# charges of the two partners. One colour per pair. Residues appear in several pairs
# (D193 in four, E194 in four), so the PAIR is coloured, not the residue.
#
#   salt      a lysine ammonium against a carboxylate
#   mgbridge  two anionic oxygens held apart by the Mg they both coordinate, NOT a salt bridge
#   hbond     a histidine or arginine nitrogen against an acceptor

dist p_K166_D193, object1 and resi 166 and name NZ, object1 and resi 193 and name OD2
dist p_K168_E194, object1 and resi 168 and name NZ, object1 and resi 194 and name OE2
dist p_K168_D193, object1 and resi 168 and name NZ, object1 and resi 193 and name OD1
dist p_K191_D193, object1 and resi 191 and name OQ2, object1 and resi 193 and name OD1
dist p_D193_E194, object1 and resi 193 and name OD1, object1 and resi 194 and name OE1
dist p_K191_E194, object1 and resi 191 and name OQ2, object1 and resi 194 and name OE1
dist p_H285_H321, object1 and resi 285 and name NE2, object1 and resi 321 and name NE2
dist p_E194_H287, object1 and resi 194 and name OE1, object1 and resi 287 and name NE2
dist p_R288_H291, object1 and resi 288 and name NE, object1 and resi 291 and name NE2
dist p_K191_H287, object1 and resi 191 and name OQ1, object1 and resi 287 and name NE2

color red, p_K166_D193
color tv_orange, p_K168_E194
color wheat, p_K168_D193
color purple, p_K191_D193
color violet, p_D193_E194
color deeppurple, p_K191_E194
color green, p_H285_H321
color limon, p_E194_H287
color forest, p_R288_H291
color palegreen, p_K191_H287

group salt, p_K166_D193 p_K168_E194 p_K168_D193
group mgbridge, p_K191_D193 p_D193_E194 p_K191_E194
group hbond, p_H285_H321 p_E194_H287 p_R288_H291 p_K191_H287

set dash_width, 3.5, salt
set dash_width, 3.5, mgbridge
set dash_width, 3.5, hbond

color deepteal, d_sc
color firebrick, d_bb
set dash_gap, 0.3
set dash_width, 2.5
set label_size, 16
set label_color, black

# ---------------------------------------------------------------- view
bg_color white
set ray_opaque_background, 1
set stick_quality, 15
set valence, 0
set ray_shadows, 0
set antialias, 2
orient ligand
zoom ligand, 6

disable object3
disable context
# The pair network is busy on top of the ligand contacts. Salt bridges are on by
# default. Type  enable mgbridge  or  enable hbond  to add the others.
disable mgbridge
disable hbond

python
from pymol import cmd
n1 = cmd.count_atoms("fix1")
n2 = cmd.count_atoms("fix2")
nl = cmd.count_atoms("ligand")
print("")
if n1 == 0 and n2 == 0 and nl == 0:
    print("  NOTHING WAS SELECTED - the script had no structure to work on.")
    print("  Objects currently loaded: %s" % cmd.get_object_list())
    print("  Load the PDB first, then run this file again:")
    print("      load /path/to/9I27_Mg_CABP.pdb")
    print("      @show_rfd3_motif.pml")
else:
    print("  object1  side-chain contributors - 11 residues - %d pinned atoms, MARINE (expect 30)" % n1)
    print("  object2  backbone contributors   -  6 residues - %d pinned atoms, SALMON (expect 13)" % n2)
    print("  object3  trans-chain chain B (E48 T53 N111) - hidden. Type: enable object3")
    print("  context  rest of the protein as faint cartoon - hidden")
    print("  ligand   CAP yellow, Mg dark sphere - %d atoms (expect 22)" % nl)
    print("")
    print("  d_sc     teal      side chain to CABP/Mg    - 11 measured contacts")
    print("  d_bb     firebrick backbone to phosphate/O4 -  5 measured contacts")
    print("")
    print("  side-chain pairs inside the motif, one colour each:")
    print("    salt      K166-D193 2.78 red | K168-E194 2.83 orange | K168-D193 3.49 wheat")
    print("    mgbridge  K191-D193 2.81 purple | D193-E194 2.88 violet | K191-E194 2.90 deeppurple")
    print("    hbond     H285-H321 2.87 green | E194-H287 3.08 limon | R288-H291 3.10 forest | K191-H287 3.48 palegreen")
    print("    salt is shown. Type  enable mgbridge  or  enable hbond  for the rest.")
    print("")
    print("  Thin grey sticks are the parts of each residue RFD3 does NOT pin.")
    if n1 != 30 or n2 != 13:
        print("")
        print("  WARNING: atom counts differ from the spec. Check chain IDs and numbering.")
print("")
python end
