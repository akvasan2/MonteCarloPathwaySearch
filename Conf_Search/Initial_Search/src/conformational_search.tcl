package require Orient
namespace import Orient::orient
namespace import ::tcl::mathfunc::min

####################################################################################################
###### Procs to create lists of theta values and phi values using fibonacci sphere with n total values. ########
####################################################################################################

proc fib_theta {n} {
	set pi 3.141592653589793
	##### gr: golden ratio in fibonacci. 
	set gr 1.618033988749895 
	
	set theta {}
	##### fill theta values according to the fibonacci relation for angles
	for {set i 0} {$i < $n} {incr i} {
		
		lappend theta [expr acos(1-2*[expr $i+1]/double($n)) ]
			
	}
	return $theta
}

proc fib_phi {n} {
	set pi 3.141592653589793
	set gr 1.618033988749895
	set phi {}
	##### fill phi values according to the fibonacci relation for angles
	for {set i 0} {$i < $n} {incr i} {
		lappend phi [expr  2*[expr $i+1]*$pi/double($gr)]
	}
	return $phi
}

############ for each theta, phi value combination, determine the xyz coordinates of the head of the axis vector ########### 

proc xyz_coor {i theta phi r} {
	set xyz {}
	set thetai [lindex $theta $i]
	set phii [lindex $phi $i]
	lappend xyz [expr $r*sin($thetai)*cos($phii)] 
	lappend xyz [expr $r*sin($thetai)*sin($phii)] 
	lappend xyz [expr $r*cos($thetai)]

	return $xyz
}

####################################################################################################
####### 1. Generating multiple antibiotic orientations using above procs ########
####################################################################################################

###################  Determine theta, phi values using fibonacci sphere of parameter num points ####################

set theta [fib_theta $num]
set phi [fib_phi $num]

################## Set a vector to use as a reference to rotate drug. #################### 
set drug [atomselect top "not protein"]

#### vector head and tail selection using head_name, tail_name parameters

set head_sel [atomselect top "not protein and name $head_name"]
set tail_sel [atomselect top "not protein and name $tail_name"]

#### vector head and tail coordinates: tv, lv 

set tv [measure center $head_sel]
set lv [measure center $tail_sel]

################## Center the middle of the vector at the origin. #################### 
set middle [vecadd $lv [vecscale 0.5 [vecsub $tv $lv]]] 
$drug moveby [vecinvert $middle]

################## Align vector to membrane normal. #################### 

set head_sel [atomselect top "not protein and name $head_name"]
set tail_sel [atomselect top "not protein and name $tail_name"]
set tv [measure center $head_sel]
set lv [measure center $tail_sel]
set A [orient $drug $tv {0 0 1}]
$drug move $A

################## Using vector generate multiple antib. orientations. #################### 

##### Use list of theta, phi values from fibonacc. sphere, as well as self rotation along vector

#### Save multi orientations to trajectory.

#### For loop over self rotations
for {set j 0} {$j<360} {incr j $self_rot} {
	$drug frame 0
	$drug update
	$drug move [transaxis z $self_rot]
	#### For loop over points in fibonacci sphere.
	for {set i 0} {$i<$num} {incr i} {
		animate dup frame 0 0
		$drug frame [expr $num*$j/double($self_rot)+$i+1]
		$drug update
		$drug move [transaxis y [lindex $theta $i] rad]
		$drug move [transaxis z [lindex $phi $i] rad]
	}
}

#### Remove initial frame (which comes from pdb structure)
animate delete beg 0 end 0 skip 0 0

#### Time to run
set showtime 1

if {$showtime==1} {
set systemTime [clock seconds]

puts "The time to generate orientations: [clock format $systemTime -format %H:%M:%S]"

}

####################################################################################################
################## 2. Translation of fib. sphere of orientations throughout protein ###################
####################################################################################################

################# Define a grid box ################################

#### set a pocket of the protein to center the grid
set pocket [atomselect top "protein and resid $pocket_residues"]
set c [measure center $pocket]

#### Define grid dimensions using del_x_pos, del_x_neg, del_y_pos, ..., del_z_neg

set max_x [expr int([lindex $c 0]) + $del_x_pos] 
set min_x [expr int([lindex $c 0]) - $del_x_neg] 

set max_y [expr int([lindex $c 1]) + $del_y_pos] 
set min_y [expr int([lindex $c 1]) - $del_y_neg] 

set max_z [expr int([lindex $c 2]) + $del_z_pos]
set min_z [expr int([lindex $c 2]) - $del_z_neg]

#### Use predefined x,y,z_spacing 

################# Generate poses 
#### Will remove any clashes with protein during this generation 

set frame 0

set nf [molinfo top get numframes]

set count [expr $nf-1]
set cclcount 0

for {set x $min_x} {$x<=$max_x} {incr x $x_spacing} {
	for {set y $min_y} {$y<=$max_y} {incr y $y_spacing} {
		for {set z $min_z} {$z<=$max_z} {incr z $z_spacing} {
			for {set f 0} {$f <$nf} {incr f} {
				incr count

				animate dup frame $f 0
				$drug frame $count
				$drug update
				$drug moveby "$x $y $z"	
				
				### Remove any poses that seriously clashes with the protein
				##  i.e. 4 clashes in which the ant-prot dist <1 Angst

				set clash [atomselect 0 "(protein within 1 of (not protein))"]
				set clash_count [$clash num]
				if {$clash_count>4} {
					animate delete beg $count end $count 0
					set count [expr $count -1]
				}
				$clash delete
			}
		}	
	}
}

animate delete beg 0 end 200 0 

##### Time to run script

set showtime 1

if {$showtime==1} {
set systemTime [clock seconds]

puts "The time is: [clock format $systemTime -format %H:%M:%S]"

}
