#########Proc to recognize if a ring piercing occurred using geometry########
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
#puts "$r"
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

proc ring_pierce_overtime {frame c5 c11 c10 c9 c8 c7 c13 c14 c15 n3 c12 c6 s1} {
# ring 1, trianle1 v0 v1 v2, triangle2 v0 v2 v3, .....
###### Atoms in any ring of the antibiotic ##############
$c5 frame $frame
$c11 frame $frame
$c10 frame $frame
$c9 frame $frame
$c8 frame $frame
$c7 frame $frame
$c13 frame $frame
$c14 frame $frame
$c15 frame $frame
$n3 frame $frame
$c12 frame $frame
$c6 frame $frame
$s1 frame $frame

####### Specify the atom content for each ring ###########

# ring 1
set v0 [measure center $c5]
set v1 [measure center $c11]
set v2 [measure center $c10]
set v3 [measure center $c9]
set v4 [measure center $c8]
set v5 [measure center $c7]


# ring 2
set v6 [measure center $c13]
set v7 [measure center $c14]
set v8 [measure center $c15]
set v9 [measure center $n3]

# ring 3
set v10 [measure center $c13]
set v11 [measure center $n3]
set v12 [measure center $c12]
set v13 [measure center $c6]
set v14 [measure center $s1]
####### Storing V0 values #########
# ring 1
dict set v0_dict 0 $v0
dict set v0_dict 1 $v0
dict set v0_dict 2 $v0
dict set v0_dict 3 $v0

# ring 2
dict set v0_dict 4 $v6
dict set v0_dict 5 $v6

# ring 3
dict set v0_dict 6 $v10
dict set v0_dict 7 $v10
dict set v0_dict 8 $v10

######### Storing V1 values ########
# ring 1
dict set v1_dict 0 $v1
dict set v1_dict 1 $v2
dict set v1_dict 2 $v3
dict set v1_dict 3 $v4

# ring 2
dict set v1_dict 4 $v7
dict set v1_dict 5 $v8

# ring 3
dict set v1_dict 6 $v11
dict set v1_dict 7 $v12
dict set v1_dict 8 $v13

######### Storing V2 values ########
# ring 1
dict set v2_dict 0 $v2
dict set v2_dict 1 $v3
dict set v2_dict 2 $v4
dict set v2_dict 3 $v5

# ring 2
dict set v2_dict 4 $v8
dict set v2_dict 5 $v9

# ring 3
dict set v2_dict 6 $v12
dict set v2_dict 7 $v13
dict set v2_dict 8 $v14


set protein [atomselect top "protein and within 2 of (not protein and name C5 C11 C10 C9 C8 C7 C13 C14 C15 N3 C12 C6 S1)" frame $frame]
set index [$protein get index]
$protein delete
set break_value 0

############ For each pose determine if a protein causes a ring piercing at that spot ################### 

foreach i $index {
        set ind_atom [atomselect top "index $i" frame $frame]
        set bonds [lindex [$ind_atom getbonds] 0]
        set bonds_close [atomselect top "index $bonds and within 2 of (not protein and name C5 C11 C10 C9 C8 C7 C13 C14 C15 N3 C12 C6 S1)" frame $frame]
        set bonds_closelist [$bonds_close get index]
        $bonds_close delete
        
        set p0 [measure center $ind_atom]
        $ind_atom delete

        foreach b $bonds_closelist {
       
                set b_atom [atomselect top "index $b" frame $frame]
                set p1 [measure center $b_atom]
                $b_atom delete

                for {set t 0} {$t<9} {incr t} {
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

#C5 C11 C10 C9 C8 C7 C13 C14 C15 N3 C12 C6 S1
set c5 [atomselect top "not protein and name C5"]
set c11 [atomselect top "not protein and name C11"]
set c10 [atomselect top "not protein and name C10"]
set c9 [atomselect top "not protein and name C9"]
set c8 [atomselect top "not protein and name C8"]
set c7 [atomselect top "not protein and name C7"]
set c13 [atomselect top "not protein and name C13"]
set c14 [atomselect top "not protein and name C14"]
set c15 [atomselect top "not protein and name C15"]
set n3 [atomselect top "not protein and name N3"]
set c12 [atomselect top "not protein and name C12"]
set c6 [atomselect top "not protein and name  C6"]
set s1 [atomselect top "not protein and name S1"]
set i 0

while {$i < [molinfo top get numframes]} {
	set out [ring_pierce_overtime $i $c5 $c11 $c10 $c9 $c8 $c7 $c13 $c14 $c15 $n3 $c12 $c6 $s1 ]
	if {$out == 1} {
		animate delete beg $i end $i 0
	} else {
		incr i		
		}
}

