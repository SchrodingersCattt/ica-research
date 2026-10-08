# ═══════════════════════════════════════════════════════════════
# NCI/IRI Auto-Render v2 — Chalk Matte + Wireframe
# Loads, sets up, exports scene, then quits VMD
# ═══════════════════════════════════════════════════════════════

if {![info exists ::env(ICA_NCI_ROOT)]} { error "Set ICA_NCI_ROOT to the local NCI input directory" }
set basedir $::env(ICA_NCI_ROOT)

# ── Display ──
display projection   Orthographic
display rendermode   GLSL
display depthcue     off
display shadows      on
display ambientocclusion on
display aoambient    0.75
display aodirect     0.30
display backgroundgradient off
color Display Background white
axes location Off
display resize 1200 1600

# Lights — softer
light 0 on
light 1 on
light 2 off
light 3 off

# ── Color scale: BGR ──
color scale method   BGR
color scale midpoint 0.666
color scale min      0.0

# Element colors
color Element Cu   orange
color Element N    blue2
color Element C    silver
color Element H    white
color Element O    red

# ── Custom materials ──
material add ChalkMatte
material change ambient   ChalkMatte 0.20
material change diffuse   ChalkMatte 0.90
material change specular  ChalkMatte 0.00
material change shininess ChalkMatte 0.00
material change mirror    ChalkMatte 0.00
material change opacity   ChalkMatte 1.00
material change outline   ChalkMatte 0.00
material change outlinewidth ChalkMatte 0.00

material add ChalkIso
material change ambient   ChalkIso 0.25
material change diffuse   ChalkIso 0.85
material change specular  ChalkIso 0.00
material change shininess ChalkIso 0.00
material change mirror    ChalkIso 0.00
material change opacity   ChalkIso 0.70
material change outline   ChalkIso 0.00
material change outlinewidth ChalkIso 0.00

# ── Load ──
mol new "$basedir/func1.cub" type cube waitfor all
mol addfile "$basedir/func2.cub" type cube waitfor all
set mol [molinfo top]
mol delrep 0 $mol

# Rep 0: atoms — chalk matte
mol representation CPK 0.40 0.15 18 16
mol color Element
mol selection {all}
mol material ChalkMatte
mol addrep $mol

# Rep 1: IRI=1.0 solid surface, colored by sign(lambda2)rho
mol representation Isosurface 1.000000 1 0 0 0 1
mol color Volume 0
mol selection {all}
mol material ChalkIso
mol addrep $mol
mol scaleminmax $mol 1 -0.04 0.02

# ── Camera ──
display resetview
rotate y by -90
rotate x by 5
scale by 1.3

# Wait for display
display update
after 3000

# ── Export scene ──
set scenefile "$basedir/scene_v2.dat"
puts "Exporting Tachyon scene to $scenefile ..."
render Tachyon "$scenefile"
puts "Scene exported. Quitting VMD."

# Quit so the batch script can continue
quit
