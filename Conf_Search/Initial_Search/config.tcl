########## Code is split into three parts: 1. Generation of multiple orientations, 2. Translation of orientations through the protein 3. Removal of ring piercings.

package require Orient
namespace import Orient::orient
namespace import ::tcl::mathfunc::min

############  psf, pdb files with protein and antibiotic ###########
##################################################################################
set psf_file "../Input_Files/system.psf"
set pdb_file "../Input_Files/system.pdb"

##################################################################################
#################### User defined parameters: ##################################
##################################################################################

############# 1. Gen of multi orientations ###############

###### number of fibonacci points on sphere used to generate theta, phi angles

set num 50 

###### head and tail atom names for vector used as reference when rotating drug

set head_name "C17" 
set tail_name "C16"

##### angle for self rotations around the vector
set self_rot 90 

############# 2. Translation throughout protein to points on a grid ###############

############# Grid creation parameters

#### pocket residues of protein the grid is centered on
set pocket_residues {16 40 42 82 132 102 106 113 114 115 116 117 118 119 120}
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

