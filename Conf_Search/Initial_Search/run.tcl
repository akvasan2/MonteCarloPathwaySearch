########## Code is split into three parts: 1. Generation of multiple orientations, 2. Translation of orientations through the protein 3. Removal of ring piercings.

package require Orient
namespace import Orient::orient
namespace import ::tcl::mathfunc::min

##################################################################################
#################### Import User defined parameters: ##################################
##################################################################################
source config.tcl

############  load psf, pdb files with protein and antibiotic ###########

mol load psf $psf_file pdb $pdb_file
puts "loaded psf, pdb"

##################################################################################
#################### Running scripts: ############################################
##################################################################################

############# 1. Gen of multi orientations and 2. Translation ###############
source conformational_search.tcl
puts "explored conformations of the drug"

############ 3. Remove ring piercings ##################################
set Ring_Names {}
for {set ring 0} {$ring < $num_rings} {incr ring} {

	lappend Ring_Names [dict get $Ring_names_dict $ring]
}

source ring_piercing.tcl

############# Save dcd files ###############################3

set frame_number [molinfo top get numframes]
set it 1

for {set t_len 0} {$t_len<$frame_number} {incr t_len 100000} {
	if {[expr $t_len+100000]>$frame_number} {
		animate write dcd dcds_test/start_${it}.dcd beg [expr $t_len+1] end -1 waitfor all
	} elseif {$tlen==0} {
		animate write dcd dcds_test/start_${it}.dcd beg [expr $t_len] end [expr $t_len+100000] waitfor all
	} else {
	animate write dcd dcds_test/start_${it}.dcd beg [expr $t_len+1] end [expr $t_len+100000] waitfor all
	}

	incr it
} 
puts "wrote dcd files"
