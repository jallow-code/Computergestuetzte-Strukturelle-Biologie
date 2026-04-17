# Orexin receptor comparison in PyMOL
# Structures:
#   6TOS chain A = OX1R with ligand NRE
#   6TPN chain A = OX2R with ligand NU8
#
# Usage:
#   pymol -r orexin_pymol.pml
# or inside PyMOL:
#   @orexin_pymol.pml
#
# After loading:
#   fig_ox1
#   fig_ox2
#   fig_overlay
#
# Render examples:
#   ray 2400,1800
#   png ox1r.png, dpi=300

reinitialize

load 6TOS.pdb, ox1_raw
load 6TPN.pdb, ox2_raw

# Receptor-only objects used for visualization
create ox1, ox1_raw and chain A and polymer.protein
create ox2, ox2_raw and chain A and polymer.protein and not resi 1000-1199

# Subtype-relevant residues for annotations
select ox1_sites, ox1 and resi 103+127
select ox2_sites, ox2 and resi 111+135
select ox1_sites_atoms, ox1_sites
select ox2_sites_atoms, ox2_sites
select pocket_focus, (ox1 and resi 103+127) or (ox2 and resi 111+135)

# Align OX2R onto OX1R using receptor backbone CA atoms
align ox2 and name CA, ox1 and name CA

# Global scene defaults for presentation-quality figures
bg_color white
set ray_opaque_background, off
set orthoscopic, on
set antialias, 2
set specular, 0.2
set shininess, 20
set ambient, 0.35
set direct, 0.55
set reflect, 0.05
set ray_shadows, on
set ray_trace_mode, 1
set ray_trace_gain, 0.08
set ray_trace_fog, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 14
set stick_radius, 0.34
set sphere_scale, 0.28
set stick_ball, on
set stick_ball_ratio, 1.6
set valence, 0
set label_font_id, 7
set label_size, 24
set label_color, black

# Color palette
set_color ox1_green, [0.24, 0.60, 0.33]
set_color ox2_gold,  [0.78, 0.60, 0.18]
set_color lig_yellow, [0.93, 0.78, 0.12]
set_color wat_blue, [0.23, 0.50, 0.88]
set_color site_cyan, [0.00, 0.62, 0.82]
set_color site_orange, [0.95, 0.45, 0.10]

color ox1_green, ox1
color ox2_gold, ox2
color site_cyan, ox1_sites
color site_orange, ox2_sites

hide everything, all
show cartoon, ox1 or ox2
show sticks, ox1_sites_atoms or ox2_sites_atoms

# Start with a clean overall orientation from the aligned receptors
orient ox1 or ox2
zoom ox1 or ox2, 8

python
from pymol import cmd

def _common():
    cmd.hide("everything", "all")
    cmd.hide("labels", "all")
    cmd.set("orthoscopic", 1)
    cmd.bg_color("white")

def fig_ox1():
    _common()
    cmd.show("cartoon", "ox1")
    cmd.show("sticks", "ox1_sites_atoms")
    cmd.show("spheres", "ox1_sites_atoms")
    cmd.color("ox1_green", "ox1")
    cmd.color("site_cyan", "ox1_sites_atoms")
    cmd.label("ox1 and name CA and resi 103", '"S103"')
    cmd.label("ox1 and name CA and resi 127", '"A127"')
    cmd.orient("ox1")
    cmd.zoom("ox1 and polymer.protein within 10 of ox1_sites", 8)

def fig_ox2():
    _common()
    cmd.show("cartoon", "ox2")
    cmd.show("sticks", "ox2_sites_atoms")
    cmd.show("spheres", "ox2_sites_atoms")
    cmd.color("ox2_gold", "ox2")
    cmd.color("site_orange", "ox2_sites_atoms")
    cmd.label("ox2 and name CA and resi 111", '"T111"')
    cmd.label("ox2 and name CA and resi 135", '"T135"')
    cmd.orient("ox2")
    cmd.zoom("ox2 and polymer.protein within 10 of ox2_sites", 8)

def fig_overlay():
    _common()
    cmd.show("cartoon", "ox1 or ox2")
    cmd.set("cartoon_transparency", 0.55, "ox1 or ox2")
    cmd.show("sticks", "ox1_sites_atoms or ox2_sites_atoms")
    cmd.show("spheres", "ox1_sites_atoms or ox2_sites_atoms")
    cmd.color("ox1_green", "ox1")
    cmd.color("ox2_gold", "ox2")
    cmd.color("site_cyan", "ox1_sites_atoms")
    cmd.color("site_orange", "ox2_sites_atoms")
    cmd.label("ox1 and name CA and resi 103", '"OX1R S103"')
    cmd.label("ox1 and name CA and resi 127", '"OX1R A127"')
    cmd.label("ox2 and name CA and resi 111", '"OX2R T111"')
    cmd.label("ox2 and name CA and resi 135", '"OX2R T135"')
    cmd.orient("pocket_focus")
    cmd.zoom("pocket_focus", 6)

cmd.extend("fig_ox1", fig_ox1)
cmd.extend("fig_ox2", fig_ox2)
cmd.extend("fig_overlay", fig_overlay)
python end

# Default to the overlay on load
fig_overlay

print("Loaded 6TOS/6TPN and aligned OX2R onto OX1R.")
print("Commands available: fig_ox1, fig_ox2, fig_overlay")
print("Render example: ray 2400,1800; png overlay.png, dpi=300")
