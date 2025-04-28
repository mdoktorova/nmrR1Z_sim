# -----------------------------------------------------
# function to return ch1 and ch2 vectors at carbon c
# of all lipids with resname $resn in the bilayer
# -----------------------------------------------------
proc vectors {resn c h1 h2} {
set result [concat]

set C [atomselect top "resname $resn and name $c"]

set HA [atomselect top "resname $resn and name $h1"]
set HB [atomselect top "resname $resn and name $h2"]

set xc [$C get x]
set yc [$C get y]
set zc [$C get z]

set xhA [$HA get x]
set yhA [$HA get y]
set zhA [$HA get z]

set xhB [$HB get x]
set yhB [$HB get y]
set zhB [$HB get z]

set xhcA [vecsub $xhA $xc]
set yhcA [vecsub $yhA $yc]
set zhcA [vecsub $zhA $zc]

set xhcB [vecsub $xhB $xc]
set yhcB [vecsub $yhB $yc]
set zhcB [vecsub $zhB $zc]

set result [concat $xhcA $yhcA $zhcA $xhcB $yhcB $zhcB]

return $result
}

# -----------------------------------------------------
# Load and analyze trajectory frame
# Frame number should be written in curr_frame.txt file
# Trajectory should be located in ../traj_center.dcd
# -----------------------------------------------------

set fnum 0
set cutoff 4

# read frame number to load
set fp [open "curr_frame.txt" r]
set frame0 [read $fp]
close $fp

# load frame
animate read dcd ../traj_center.dcd beg $frame0 end $frame0 waitfor all

# -----------------------------------------------------
# read information about the CH bonds
# file ch_vectors.txt should be organized as follows:
# 1st line: lipid resname (e.g. DMPC)
# 2nd, 3rd, etc. line: [carbon name] [hydrogen1 name] [hydrogen2 name]
#
# *NOTE* for double bonds (i.e. carbon atom with only 1 hydrogen),
# write the name of the hydrogen atom twice
# -----------------------------------------------------

set fp [open "ch_vectors.txt" r]
set data [read $fp]
close $fp

set resn [lindex $data 0]
set nvecs [expr ([llength $data]-1)/3.0]

# -----------------------------------------------------
# go through all CH vectors (equal in number to the number of lines
# in ch_vectors.txt minus 1)
# for each CH vector, get its direction and append it to
# the appropriate carbon-specific file
# -----------------------------------------------------
for {set i 1} {$i<=$nvecs} {incr i} {
    set c [lindex $data [expr ($i-1)*3+1]]
    set h1 [lindex $data [expr ($i-1)*3+2]]
    set h2 [lindex $data [expr ($i-1)*3+3]]

    set res [vectors $resn $c $h1 $h2]
    set outf [open angle_$c.txt a+]
    puts $outf $res
    close $outf
}

exit

