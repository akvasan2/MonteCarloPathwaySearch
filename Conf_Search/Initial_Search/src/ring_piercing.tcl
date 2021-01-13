#########Proc to recognize if a ring piercing occurred using geometrical analysis########
# Code modified from http://geomalgorithms.com/a06-_intersect-2.html 
# Copyright 2001 softSurfer, 2012 Dan Sunday
#############################################################################

proc intersect3d_raytriangle { p0 p1 v0 v1 v2 } {

set u [vecsub $v1 $v0]
set v [vecsub $v2 $v0]
set n [veccross $u $v]

set dir [vecsub $p1 $p0]
set w0 [vecsub $p0 $v0]
set a [expr {-1 * [vecdot $n $w0] } ]
set b [vecdot $n $dir]

if {abs($b)<0.01} {
	return 0		
}

set r [expr {$a/double($b)}]

if {$r<0} {
	return 0
}

if {$r>1.0} {
	return 0
}

set I_st [vecadd $p0 [vecscale $r $dir]]

set uu [vecdot $u $u]
set uv [vecdot $u $v]
set vv [vecdot $v $v]
set w [vecsub $I_st $v0]
set wu [vecdot $w $u]
set wv [vecdot $w $v]
set D [expr {$uv*$uv - $uu*$vv}]
set s [expr {($uv*$wv - $wu*$vv)/double ($D)}]

if {$s <0 || $s>1} {

	return 0
}

set t [expr {($uv*$wu - $uu*$wv)/double ($D)}]

if {$t<0 || [expr $s +$t ] >1 } {
	return 0
}


return 1

}

########Proc to run ring piercing script over pose data ############
#####################################################################

proc ring_pierce_overtime {frame Ring_Names Ring_sel} {
# ring 1, triangle1 v0 v1 v2, triangle2 v0 v2 v3, .....
###### Atoms in any ring of the antibiotic ##############

set Ring_Names_unique_prelim {}

for {set ring 0 } {$ring<[llength $Ring_Names]} {incr ring} {
	set sels {}
	for {set r_ind 0} {$r_ind < [llength [lindex $Ring_Names $ring]]} {incr r_ind } {
		lappend Ring_Names_unique_prelim [lindex [lindex $Ring_Names $ring] $r_ind]
	}

}

set Ring_Names_unique [lsort -unique $Ring_Names_unique_prelim]



for {set ring 0} {$ring <[llength $Ring_sel]} {incr ring} { 
	for {set r_ind 0} {$r_ind < [llength [lindex $Ring_sel $ring]]} {incr r_ind } {

		[lindex [lindex $Ring_sel $ring] $r_ind] frame $frame  

		}	
	}	

set Ring_at_cents {}
for {set ring 0} {$ring <[llength $Ring_sel]} {incr ring} { 
	set cent_list {}
	for {set r_ind 0} {$r_ind < [llength [lindex $Ring_sel $ring]]} {incr r_ind } {

		lappend cent_list [measure center [lindex [lindex $Ring_sel $ring] $r_ind]]  
		#puts $cent_list
		}	
		lappend Ring_at_cents $cent_list
	}	


####### Create dictionaries storing vertex information #########

####### Storing V0 values #########

set v0_it -1

for {set ring 0} {$ring <[llength $Ring_sel]} {incr ring} { 
	
	for {set r_ind 0} {$r_ind < [expr [llength [lindex $Ring_sel $ring]] - 2]} {incr r_ind } {
		incr v0_it	

		dict set v0_dict $v0_it [lindex [lindex $Ring_at_cents $ring] 0]  
		}

		
	}	

####### Storing V1 values #########

set v1_it -1

for {set ring 0} {$ring <[llength $Ring_sel]} {incr ring} { 
	
	for {set r_ind 1} {$r_ind < [expr [llength [lindex $Ring_sel $ring]] - 1]} {incr r_ind } {
		incr v1_it		

		dict set v1_dict $v1_it [lindex [lindex $Ring_at_cents $ring] $r_ind]  
		}

		
	}	

####### Storing V2 values #########
set v2_it -1

for {set ring 0} {$ring <[llength $Ring_sel]} {incr ring} { 
	
	for {set r_ind 2} {$r_ind < [expr [llength [lindex $Ring_sel $ring]]]} {incr r_ind } {
		incr v2_it	

		dict set v2_dict $v2_it [lindex [lindex $Ring_at_cents $ring] $r_ind]  
		}

		
	}	


set protein [atomselect top "protein and within 2 of (not protein and name $Ring_Names_unique )" frame $frame]
set index [$protein get index]
$protein delete
set break_value 0

############ For each pose run script to determine if the protein causes a ring piercing at that spot ################### 

foreach i $index {
        set ind_atom [atomselect top "index $i" frame $frame]
        set bonds [lindex [$ind_atom getbonds] 0]
        set bonds_close [atomselect top "index $bonds and within 2 of (not protein and name $Ring_Names_unique )" frame $frame]
        set bonds_closelist [$bonds_close get index]
        $bonds_close delete
        
        set p0 [measure center $ind_atom]
        $ind_atom delete

        foreach b $bonds_closelist {
       
                set b_atom [atomselect top "index $b" frame $frame]
                set p1 [measure center $b_atom]
                $b_atom delete
                ### iterate over all triangles (can be described by number of $v2_it)
                for {set t 0} {$t <= $v2_it} {incr t} {
                        set v0_ind [dict get $v0_dict $t]
                        set v1_ind [dict get $v1_dict $t]
                        set v2_ind [dict get $v2_dict $t]
                        set intersect [intersect3d_raytriangle $p0 $p1 $v0_ind $v1_ind $v2_ind]
                        if {$intersect == 1} {
                                puts "Broke triangles in frame $frame and tri $t"
                                return 1
                        }

                }
        }
}
return 0
}


############ Running the ring piercing script #############

set Ring_sel {}
for {set ring 0 } {$ring<[llength $Ring_Names]} {incr ring} {
	set sels {}
	for {set r_ind 0} {$r_ind < [llength [lindex $Ring_Names $ring]]} {incr r_ind } {
		set at_name [lindex [lindex $Ring_Names $ring] $r_ind]
		lappend sels [atomselect top "not protein and name $at_name"]
	}
lappend Ring_sel $sels

}

set i 0

while {$i < [molinfo top get numframes]} {
	set out [ring_pierce_overtime $i $Ring_Names $Ring_sel]
	if {$out == 1} {
		animate delete beg $i end $i 0
		
	} else {
		incr i		
		}
}

