package require Orient
namespace import Orient::orient
namespace import ::tcl::mathfunc::min

##################################################################################
#################### User defined parameters: ##################################
##################################################################################

############  psf, pdb files with protein and antibiotic ###########
set psf_file "../Input_Files/system.psf"
set pdb_file "../Input_Files/system.pdb"

############# 1. Gen of multi orientations ###############

###### number of fibonacci points on sphere used to generate theta, phi angles
set num 2 

###### head and tail atom names for vector used as reference when rotating drug
set head_name "C17" 
set tail_name "C16"

##### angle for self rotations around the vector
set self_rot 90 

############# 2. Translation throughout protein to points on a grid ###############

############# Grid creation parameters

#### pocket residues of protein the grid is centered on
set pocket_residues {16 40 42 82 132 102 106 113 114 115 116 117 118 119 120}
#
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

############# 3. Ring piercings ###############

### Number of  rings in system
set num_rings 3

### manually set list of atom names for each ring in system

dict set Ring_names_dict 0 {C1 C2 C3 C4 C5 C6}
dict set Ring_names_dict 1 {C11 C12 C13 S1 N2}
dict set Ring_names_dict 2 {C9 C10 C11 N2}

##################################################################################
#################### Running scripts: ############################################
##################################################################################

############  load psf, pdb files with protein and antibiotic ###########
mol load psf $psf_file pdb $pdb_file
puts "loaded psf, pdb"

############# 1. Gen of multi orientations and 2. Translation ###############
source src/conformational_search.tcl
puts "explored conformations of the drug"

############ 3. Remove ring piercings ##################################
set Ring_Names {}
for {set ring 0} {$ring < $num_rings} {incr ring} {

	lappend Ring_Names [dict get $Ring_names_dict $ring]
}

source src/ring_piercing.tcl

############# Save dcd files ###############################3

set frame_number [molinfo top get numframes]
set it 1

for {set t_len 0} {$t_len<$frame_number} {incr t_len 100000} {
	if {[expr $t_len+100000]>$frame_number} {
		animate write dcd dcds/dcdsstart_${it}.dcd beg [expr $t_len+1] end -1 waitfor all
	} elseif {$t_len==0} {
		animate write dcd dcds/dcdsstart_${it}.dcd beg [expr $t_len] end [expr $t_len+100000] waitfor all
	} else {
	animate write dcd dcds/dcdsstart_${it}.dcd beg [expr $t_len+1] end [expr $t_len+100000] waitfor all
	}

	incr it
} 
puts "wrote dcd files"
