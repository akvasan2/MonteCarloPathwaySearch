###### Script to determine z, inclination, and azimuthal angles, which are slow coordiantes for antibiotic permeation#######
#########################################################################
##################Important parameters to set: ##################
# head_name and tail_name: head and tail of reference vector for calculation.
#  output_file: output file name
#  num_dcds: number of minimized pose data dcds

proc drug_vec { frame head_name tail_name } {
    set head [atomselect top "not protein and name $head_name" frame $frame]
    set tail [atomselect top "not protein and name $tail_name" frame $frame]
    set tv [measure center $head]
    set lv [measure center $tail]
    return [vecnorm [vecsub $tv $lv]]
}


##### parameters to run calculations ######

### vector head and tail
set head_name "C17"
set tail_name "C16"
set output_file "slow_coor.dat"
### total number of pose dcds
set num_dcds 8

set file [open $output_file w]

##### Running calculation ############ 

for {set j 1} {$j <= $num_dcds} {incr j} {
mol load psf ../../Input_Files/system.psf
mol addfile ../../Minimization/minimize_output/dcdsstart_$j.dcd waitfor all

set sel_all [atomselect top "not protein"]

set nf [molinfo top get numframes]

for {set i 0} {$i < $nf} {incr i} {
	set drug_vec [ drug_vec $i $head_name $tail_name ]
	#### evaluate theta and phi angles using geometrical relations
	set theta [expr {acos ([ lindex $drug_vec 2])}]
	### Since phi can be between 0 and 2pi in radians need to apply this "trick"

	if {[lindex $drug_vec 0] < 0} {
	set phi [expr {atan ([lindex $drug_vec 1]/[lindex $drug_vec 0]) + 3.14159265}]
	puts $phi
	} else {
	set phi [expr {atan ([lindex $drug_vec 1]/[lindex $drug_vec 0])}]
	}
	if {$phi < 0} {
	set phi [expr $phi + 2*3.14159265]
	}

	$sel_all frame $i
	puts $file "[lindex [measure center $sel_all] 2] [expr $theta*180/3.14159265] [expr $phi*180/3.14159265]"
	}
mol delete all
}

close $file

