proc drug_vec { frame } {
    set C22 [atomselect top "not protein and name C17" frame $frame]
    set C7 [atomselect top "not protein and name C16" frame $frame]
    set tv [measure center $C22]
    set lv [measure center $C7]
    return [vecnorm [vecsub $tv $lv]]
}

set file [open "angles_z.dat" w]
for {set j 1} {$j < 9} {incr j} {
mol load psf ../../creating_dcds/system.psf
mol addfile ../../minimize_GBIS_small_XY/minimize_output/Min_200/dcdsstart_$j.dcd waitfor all

set sel_all [atomselect top "not protein"]

set nf [molinfo top get numframes]


for {set i 0} {$i < $nf} {incr i} {
	set drug_vec [drug_vec $i]
	set theta [expr {acos ([ lindex $drug_vec 2])}]
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

