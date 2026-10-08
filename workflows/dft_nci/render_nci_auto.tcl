# ═══════════════════════════════════════════════
# Auto-render NCI/IRI isosurface with VMD + Tachyon
# Run with: vmd.exe -e render_nci_auto.tcl
# (Must be in GUI mode, not text mode)
# ═══════════════════════════════════════════════

if {![info exists ::env(ICA_NCI_ROOT)]} { error "Set ICA_NCI_ROOT to the local NCI input directory" }
set basedir $::env(ICA_NCI_ROOT)
set tachyon ""
if {[info exists ::env(ICA_TACHYON_EXE)]} { set tachyon $::env(ICA_TACHYON_EXE) }

# Display
display projection   Orthographic
display rendermode   GLSL
display depthcue     off
display shadows      on
display ambientocclusion on
display aoambient    0.80
display aodirect     0.25
display backgroundgradient off
color Display Background white
axes location Off
display resize 1200 1600

# Lights
light 0 on
light 1 on
light 2 on
light 3 on

# Color scale
color scale method   BGR
color scale midpoint 0.666
color scale min      0.0

# Element colors
color Element Cu orange
color Element N  blue
color Element C  gray
color Element H  white

# Load
mol new "$basedir/func1.cub" type cube waitfor all
mol addfile "$basedir/func2.cub" type cube waitfor all
set mol [molinfo top]
mol delrep 0 $mol

# Rep 0: Small ball-and-stick
mol representation CPK 0.40 0.15 18 16
mol color Element
mol selection {all}
mol material Diffuse
mol addrep $mol

# Rep 1: IRI isosurface at 1.0, colored by sign(λ₂)ρ
mol representation Isosurface 1.000000 1 0 0 1 1
mol color Volume 0
mol selection {all}
mol material EdgyShiny
mol addrep $mol
mol scaleminmax $mol 1 -0.04 0.02

# Material tweaks
material change ambient   Diffuse  0.15
material change diffuse   Diffuse  0.85
material change specular  Diffuse  0.05
material change opacity   Diffuse  1.00

material change ambient   EdgyShiny 0.10
material change diffuse   EdgyShiny 0.55
material change specular  EdgyShiny 0.60
material change shininess EdgyShiny 0.70
material change opacity   EdgyShiny 0.82

# Camera: look along A-axis
display resetview
rotate y by -90
rotate x by 5
scale by 1.3

# Wait for display to settle
after 2000

# Auto render with Tachyon
set scenefile "$basedir/scene.dat"
set outfile   "$basedir/NCI_3d_vmd.tga"

puts "Rendering scene to $scenefile ..."
render Tachyon "$scenefile"

puts "Running Tachyon ray tracer..."
set cmd [list "$tachyon" "$scenefile" -aasamples 12 -res 4000 5200 -format TGA -o "$outfile"]
puts "CMD: $cmd"

if {[catch {exec {*}$cmd} result]} {
    puts "Tachyon error: $result"
    puts "Trying without quotes..."
    catch {exec "$tachyon" "$scenefile" -aasamples 12 -res 4000 5200 -format TGA -o "$outfile"} result2
    puts "Result: $result2"
}

puts ""
puts "═══════════════════════════════════"
puts "DONE! Output: $outfile"
puts "═══════════════════════════════════"

# Don't quit — let user inspect in VMD
# quit
