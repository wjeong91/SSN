# RFD3 motif for the 9I27 CABP site, split by how each residue reaches the ligand.
#
#   object1  side-chain contributors  K166 K168 KCX191 D193 E194 H285 H287 R288 H291 H321 K329
#   object2  backbone contributors    S368 G369 G370 A392 G393 G394
#   object3  trans-chain, chain B     E48 T53 N111   (only in the 7-component spec; starts hidden)
#
# Within object1 and object2 the atoms RFD3 actually pins are drawn thick and in colour;
# the rest of each residue is thin and grey. Site shown is CAP A501 / MG A500, built from
# the chain A C-terminal domain plus the chain B N-terminal domain.
# Every distance drawn is one measured in RFD3_ATOM_MEASUREMENTS.md.
#
# run:  pymol 9I27_Mg_CABP.pdb show_rfd3_motif.pml
#   or, inside PyMOL:  @show_rfd3_motif.pml

python
from pymol import cmd
if "9I27_Mg_CABP" not in cmd.get_object_list():
    cmd.load("9I27_Mg_CABP.pdb", "9I27_Mg_CABP")
python end

hide everything
remove solvent

# ---------------------------------------------------------------- objects
create object1, 9I27_Mg_CABP and chain A and resi 166+168+191+193+194+285+287+288+291+321+329
create object2, 9I27_Mg_CABP and chain A and resi 368+369+370+392+393+394
create object3, 9I27_Mg_CABP and chain B and resi 48+53+111
create ligand,  9I27_Mg_CABP and chain A and resn CAP+MG
create context, 9I27_Mg_CABP and polymer

# The 30 + 13 atoms the spec pins.
select fix1, object1 and ( \
      (resi 166+168+329 and name CE+NZ) \
   or (resi 191 and name CX+OQ1+OQ2) \
   or (resi 193 and name CG+OD1+OD2) \
   or (resi 194 and name CD+OE1+OE2) \
   or (resi 285+287+291+321 and name ND1+CE1+NE2) \
   or (resi 288 and name NE+CZ+NH2))
select fix2, object2 and ( \
      (resi 368 and name C+O) \
   or (resi 369+393 and name N+CA+C) \
   or (resi 370+394 and name N+CA) \
   or (resi 392 and name C))

# ---------------------------------------------------------------- style
show sticks, object1 or object2 or object3
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
color magenta, ligand and resn MG

# Category colours go on the pinned atoms only.
# Category colour covers the whole pinned atom set, heteroatoms included, so the
# side-chain / backbone split is readable at a glance. For CPK heteroatoms instead,
# run:  util.cnc("fix1 or fix2")
color marine, fix1
color salmon, fix2
color palegreen, object3

# ---------------------------------------------------------------- measured contacts
# Side chain to CABP or Mg.
dist d_sc, object1 and resi 166 and name NZ,  ligand and name O2
dist d_sc, object1 and resi 168 and name NZ,  ligand and name O6
dist d_sc, object1 and resi 191 and name OQ2, ligand and resn MG
dist d_sc, object1 and resi 191 and name OQ1, ligand and name O3
dist d_sc, object1 and resi 193 and name OD1, ligand and resn MG
dist d_sc, object1 and resi 194 and name OE1, ligand and resn MG
dist d_sc, object1 and resi 287 and name NE2, ligand and name O3
dist d_sc, object1 and resi 288 and name NE,  ligand and name O5P
dist d_sc, object1 and resi 288 and name NH2, ligand and name O4P
dist d_sc, object1 and resi 321 and name ND1, ligand and name O6P
dist d_sc, object1 and resi 329 and name NZ,  ligand and name O2P

# Backbone to the phosphates and O4.
dist d_bb, object2 and resi 368 and name O,  ligand and name O6P
dist d_bb, object2 and resi 369 and name N,  ligand and name O4
dist d_bb, object2 and resi 370 and name N,  ligand and name O2P
dist d_bb, object2 and resi 393 and name N,  ligand and name O3P
dist d_bb, object2 and resi 394 and name N,  ligand and name O1P

# The two second-shell braces, both 1.00 across every CABP-bound form II structure.
dist d_brace, object1 and resi 285 and name NE2, object1 and resi 321 and name NE2
dist d_brace, object1 and resi 291 and name NE2, object1 and resi 288 and name NE

color deepteal,  d_sc
color firebrick, d_bb
color orange,    d_brace
set dash_gap, 0.3
set dash_width, 2.5
set label_size, 16
set label_color, black
set label_position, (0, 0, 1.5)

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

python
print("")
print("  object1  side-chain contributors - 11 residues - 30 pinned atoms in MARINE")
print("  object2  backbone contributors   -  6 residues - 13 pinned atoms in SALMON")
print("  object3  trans-chain chain B (E48 T53 N111) - hidden. Type: enable object3")
print("  context  rest of the protein as faint cartoon - hidden")
print("")
print("  d_sc     teal      side chain to CABP/Mg    - 11 measured contacts")
print("  d_bb     firebrick backbone to phosphate/O4 -  5 measured contacts")
print("  d_brace  orange    H285-H321 2.87 and H291-R288 3.10")
print("")
print("  Thin grey sticks are the parts of each residue RFD3 does NOT pin.")
print("")
python end
