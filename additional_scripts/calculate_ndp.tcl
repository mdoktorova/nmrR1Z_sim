package require density_profile
set tcl_precision 5

# load trajectory to be analyzed
animate read dcd ../traj_center.dcd

# get the total number of frames
set nf [molinfo top get numframes]
set lastFrame [expr $nf-1]

set frameFROM 0
set frameTO $lastFrame

# set the slabs nmber/dimensions (z-bins)
set dp [density_profile -rho number -selection "all" -frame_to $lastFrame -frame_from $lastFrame -average 1 -resolution 0.2]
set dp_f [list $dp]
set bin_names [lindex $dp 1]
set nbins [llength [lindex $dp 0]]

# total number of NDPs to be calculated
set atomnfull 1

# get NDP of while bilayer; edit selection to contain the resnames of all lipids in the bilayer
set dp [density_profile -rho number -selection "resname DLPC CHL1" -frame_to $frameTO -frame_from $frameFROM -average 1 -resolution 0.2]
set dp_f [concat $dp_f [list $dp]]

# write the profile to file
set file [open "ndp.dat" w]
for {set bin 0} {$bin < $nbins} {incr bin} {
    set row [lindex $bin_names $bin]
    for {set atom 1} {$atom <= $atomnfull} {incr atom} {
	if {[lindex $dp_f $atom 1 0] > [lindex $bin_names $bin] ||
	    [lindex $dp_f $atom 1 [expr [llength [lindex $dp_f $atom 1]]-1]] < [lindex $bin_names $bin]} {
	    set row [concat $row 0.0]
	    continue
	}

	for {set ind 0} {$ind < [llength [lindex $dp_f $atom 1]]} {incr ind} {
	    if {[lindex $dp_f $atom 1 $ind] == [lindex $bin_names $bin]} {
		set row [concat $row [lindex $dp_f $atom 0 $ind]]
		break
	    }
	}
    }
    puts $file $row
}

close $file

exit



