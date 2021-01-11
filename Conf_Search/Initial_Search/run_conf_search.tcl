########## Code is split into three parts: 1. Generation of multiple orientations, 2. Translation of orientations through the protein 3. Removal of ring piercings.

package require Orient
namespace import Orient::orient
namespace import ::tcl::mathfunc::min

############  load psf, pdb files with protein and antibiotic ###########
##################################################################################
set psf_file "system.psf"
set pdb_file "system.pdb"
mol load psf $psf_file pdb $pdb_file
puts "loaded psf, pdb"

#################### Definition of parameters: ##################################
##################################################################################

############# 1. Gen of multi orientations ###############
##################################################################################

###### number of fibonacci points on sphere used to generate theta, phi angles

###### number of fibonacci points
set num 50 

###### head and tail atom names for vector used as reference when rotating drug

set head_name "C17" 
set tail_name "C16"

##### angle for self rotations around the vector
set self_rot 90 

############# 2. Translation throughout protein to points on a grid ###############
##################################################################################

############# Grid creation parameters

#### pocket of protein the grid is centered on
set pocket [atomselect top "protein and resid 16 40 42 82 132 102 106 113 114 115 116 117 118 119 120"] 
#### Grid dimensions: positive and negative max distance from center

set del_x_pos 8
set del_x_neg 8

set del_y_pos 8
set del_y_neg 8

set del_z_pos 20
set del_z_neg 14

#### Grid spacing in each dimension
set x_spacing 1
set y_spacing 1
set z_spacing 1

##################################################################################

##################################################################################
#################### Running scripts: ############################################
##################################################################################

############# 1. Gen of multi orientations and 2. Translation ###############
source conformational_search.tcl
puts "explored conformations of the drug"

############ 3. Removal of ring piercings ##################################
#source ring_piercing.tcl

############# Save dcd files ###############################3
set frame_number [molinfo top get numframes]
set it 1

for {set t_len 0} {$t_len<$frame_number} {incr t_len 100000} {
	if {[expr $t_len+100000]>$frame_number} {
		animate write dcd dcds/start_${it}.dcd beg [expr $t_len+1] end -1 waitfor all
	} elseif {$tlen==0} {
		animate write dcd dcds/start_${it}.dcd beg [expr $t_len] end [expr $t_len+100000] waitfor all
	}
	else {
	animate write dcd dcds/start_${it}.dcd beg [expr $t_len+1] end [expr $t_len+100000] waitfor all
	}

	incr it
} 
puts "wrote dcd files"
