###### File used to make fix.pdb and PIE.pdb files used for minimization and pair interaction energy calculation

mol load psf system.psf pdb system.pdb

#### making PIE.pdb

set sel1 [atomselect top "protein"]
set sel2 [atomselect top "not protein"]
set all [atomselect top all]
$all set beta 0
$sel1 set beta 1
$sel2 set beta 2
$all writepdb PIE.pdb

#### making fix.pdb

$all set occupancy 0
$all set beta 0
set sel [atomselect top "protein and backbone"]
set all [atomselect top all]
$all set occupancy 0
$sel set occupancy 1
$all writepdb fix.pdb
