####################### Using the path_frames.dat file, choose the best possible snapshot at each step of the path #################

## proc to measure distance between two objects projected in x-y plane
proc lat_dist { xval_old yval_old xval_new yval_new } {
        return [expr { (($xval_new - $xval_old)**2 + ($yval_new-$yval_old)**2)**0.5  } ] 
}

## proc to obtain a random index from an array

proc lrandom {L} {
    return [lindex $L [expr {int(rand()*[llength $L])}]]
}

######## User defined parameters ###########

## number of dcds stored
set num_dcds 8

###Cutoff values for how much the transition, rotational, and conformational (internal degrees of freedom) deviation can be in each step of the algorithm

set trans_cutoff 5
set rot_cutoff 5
set conf_cutoff 5

mol new ../../../../Conf_Search/Input_Files/system.psf

for {set dcd_it 1} {$dcd_it <= num_dcds} {incr $dcd_it} {

	mol addfile ../../../Conf_Search/Minimization/minimize_output/dcdsstart_${dcd_it}.dcd
}

set data_id [molinfo top get id]

puts "Loaded data"

## File containing all possible frames at each step of the path
set fil [open path_frames.dat]
set lines [split [read $fil] "\n"]
close $fil

## Number of steps in the path:
set numlines [expr [llength $lines] -1 ]
puts $numlines

## file to save the final accepted path
set path_file [open accepted_pathway.dat "w"]

################################################################
########## Running loop to obtain an appropriate path.##########
################################################################

set frames_top [lindex $lines 0]
set accepted 0

### We try to get an appropriate path by running a while loop until an accepted path is reached ###
while {$accepted ==0} {
	## initialize list of accepted frames in the path
        set accepted_frmes {}

        ## start at a random initial frame at the beginning of the path.
	set iniframe [lrandom $frames_top ]
	set iniframe_sel [atomselect $data_id "all" frame $iniframe]
	$iniframe_sel writepdb PDBs/window_0.pdb
	$iniframe_sel delete
        
	lappend accepted_frmes $iniframe
        
	## Now go to each step of the path to choose an appropriate frame
        for  {set l 1} {$l<[expr $numlines + 1 ]} {incr l} {        
		
		## fo: accepted frame before this particular step
                set fo [lindex $accepted_frmes [expr $l -1 ]]
                set fo [expr int($fo)]
		
		## atomselection of drug at this previous step
		mol new PDBs/window_[expr $l -1 ].pdb
		set old_id [molinfo top get id]
                set drug_old [atomselect $old_id "not protein"]
		
		#frame $fo]

		## center of drug at previous step
                set centold [measure center $drug_old]
                set xold [lindex $centold 0]
                set yold [lindex $centold 1]
                
		## initialize array to contain all frames within translational cutoff
                set possible_frmes_t {}
		## initialize array to contain all frames within rotational cutoff
                set possible_frmes_r {}
                set min_rmsd_conf 30000
                set possible_frme_c 0
                set linenew [lindex $lines $l]
                ## loop over list of possible frames at this step of the path to see if it is translationally connected
                foreach fr $linenew {
                        set fr_i [expr int($fr)]
                 
                        set drug_new [atomselect $data_id "not protein" frame $fr_i]
                 
			## measure x and y coordinates of this frame
                        set centnew [measure center $drug_new]
                        set xnew [lindex $centnew 0]
                        set ynew [lindex $centnew 1]
                 	
			## only accept to list of allowed frames if x and y are within cutoff dist to previous step
                        set dist [lat_dist $xold $yold $xnew $ynew] 
                        if {$dist < $trans_cutoff} {
                                lappend possible_frmes_t $fr_i
                        }
                        $drug_new delete
                }        
		## if no possible frames were determined here then break loop and try to find a new path
                if {[llength $possible_frmes_t]==0} {
                       puts "broken translation"
                       mol delete $old_id
                       break 
                }
		#puts "Translational loop done"
		## loop over list of translational accepted frames to see if they are rotationally connected to previous step

                foreach fs $possible_frmes_t {
                        set fs_i [expr int($fs)]
                        set drug_new [atomselect $data_id "not protein" frame $fs_i]
			## align center of previous step drug to trial frame. measure RMSD to obtain estimate for change in rotation 
                        set centnew [measure center $drug_new]
                        set dist_newold [vecsub $centnew $centold]
                        $drug_old moveby $dist_newold
                        set rmsdval [measure rmsd $drug_old $drug_new]
			## move old drug back to original location
                        $drug_old moveby [vecinvert $dist_newold]
			## if rmsd < cutoff then add to list of acceptable frames connected in rotational space 
                        if { $rmsdval < $rot_cutoff } {        
                                lappend possible_frmes_r $fs_i 
                        }
                        $drug_new delete
                }        
		## if the length of acceptable frames here is 0 then go to start and try to find new path

                if {[llength $possible_frmes_r]==0} {
                       puts "broken rotation"
                       mol delete $old_id
                       puts $l
                       break 
                }
		## using list of filtered poses within transl, rotat cutoffs, find poses with small conformational deviation from prev. step
		#puts "Rotational loop done"
                foreach fc $possible_frmes_r {
                        set fc_i [expr int($fc)]
                        set drug_new [atomselect $data_id "not protein" frame $fc_i]
			## align frame from previous step to trial frame
                        set trans_matrix [measure fit $drug_old $drug_new]
                        $drug_old move $trans_matrix
			## measure rmsd difference between trial frame and previous pose 
                        set rmsdval [measure rmsd $drug_old $drug_new]
			## if rmsd is less than predefined min_rmsd_conf value then replace min_rmsd_conf and possible_frme_c
                        if {$rmsdval < $min_rmsd_conf} {
                                set min_rmsd_conf $rmsdval
                                set possible_frme_c $fc_i
                        }
                        $drug_new delete 
                 }
		 set drug_step [atomselect $data_id "all" frame $possible_frme_c ]
		 $drug_step writepdb PDBs/window_$l.pdb
		 $drug_step delete
		 puts "Obtained ${l} step in path "
		 mol delete $old_id
		 ## at the end of the loop, $possible_frme_c is the accepted frame and you continue with the next step in the path
                 lappend accepted_frmes $possible_frme_c 
        }        
	## if you find a path that has all frames, then this is your final path! the while loop is then finished
        if {[llength $accepted_frmes ] >= 45 } { 
               puts "length of frames is" 
               puts [llength $accepted_frmes] 
               puts $accepted_frmes
               puts $path_file $accepted_frmes 
               incr accepted 
               puts "accepted"        
        }
}

close $path_file 

