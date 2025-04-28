# count number of frames
set nf [molinfo top get numframes]

# open file for writing
set file [open box2.txt w]

# initialize variables to store the mean values
# of the box dimensions
set meanX 0
set meanY 0
set meanZ 0

# loop through the frames
for {set i 0} {$i<$nf} {incr i} {

    # go to frame i
    animate goto $i

    # get (x,y,z)=(a,b,c) dimensions of the box
    set a [molinfo top get a]
    set b [molinfo top get b]
    set c [molinfo top get c]

    # write them to the file, together with the frame number
    puts $file "$i $a $b $c"

    # update meanX, meanY and meanZ
    set meanX [expr $meanX+(1.0/$nf)*$a]
    set meanY [expr $meanY+(1.0/$nf)*$b]
    set meanZ [expr $meanZ+(1.0/$nf)*$c]
}

close $file

# output the mean box dimensions
puts [format "mean (x,y,z) dimensions: (%.2f, %.2f, %.2f)" $meanX $meanY $meanZ]
