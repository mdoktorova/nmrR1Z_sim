package require pbctools

# get total number of frames
set nf [molinfo top get numframes]

# first frame is 0, so last frame is nf-1
set nff [expr $nf-1]

# center the trajectory on one terminal methyl group of one of the lipids
# this is necessary to avoid potential problems with the actual centering later
pbc wrap -centersel "name C218 and resid 1" -center com -compound residue -all

# set the center of the bilayer midplane to be at (x,y,z)=(0,0,0)
for {set i 0} {$i<$nf} {incr i} {

    # go to frame $i
    animate goto $i

    # calculate the geometric center of all methyl groups
    # modify selection appropriately to contain the methyl groups of all lipids
    set terminal_methyls [atomselect top "resname POPC and name C218 C316" frame $i]
    set ctr [measure center $terminal_methyls]

    # select all atoms in the system
    set all [atomselect top "all" frame $i]

    # move all atoms in x, y and z so that the center of the methyl groups
    # is at (0,0,0)
    set offset [expr 0-[lindex $ctr 0]]
    set off [list $offset 0 0]
    $all moveby $off

    set offset [expr 0-[lindex $ctr 1]]
    set off [list 0 $offset 0]
    $all moveby $off
    
    set offset [expr 0-[lindex $ctr 2]]
    set off [list 0 0 $offset]
    $all moveby $off
}

# center the whole box on (x,y,z)=(0,0,0)
pbc wrap -center origin -compound residue -all

# write the centered trajectory to a new file
animate write dcd traj_center.dcd waitfor all

exit 
