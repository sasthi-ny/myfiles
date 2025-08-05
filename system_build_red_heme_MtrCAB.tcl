package require psfgen

resetpsf

mol new final_mem_mtrAB_SO_withca_ions.psf
mol addfile equilibration_1_5.coor
set sel1 [atomselect top "segname PA PB HEMA"]
set sel2 [atomselect top "all not segname PA PB HEMA"]
$sel1 writepdb protein_heme.pdb
$sel1 writepsf protein_heme.psf

$sel2 writepdb mem-wat-ion-ca.pdb
$sel2 writepsf mem-wat-ion-ca.psf
$sel1 delete
$sel2 delete

mol delete all

## Loading topology files
topology ../../../toppar/top_all36_prot.rtf
topology ../../../toppar/toppar_water_ions.str
topology ../../../toppar/top_hemeRed_charmm36.inp
## End loading topology files

mol new protein_heme.psf
mol addfile protein_heme.pdb

#pdbalias residue HEC HER

set all [atomselect top all]

set cname [atomselect top "segname PA and resid 67 70 100 103 136 139 160 163 183 186 209 212 233 236 255 258 278 281 313 316"]
set hname [atomselect top "segname PA and resid 71 85 104 110 140 153 164 167 187 200 213 216 237 248 259 262 282 287 317 321"]
set heme [atomselect top "resname HEC"]
$cname set resname CYR
$hname set resname HSR
$heme set resname HER
$all writepdb protein_heme_mod.pdb
#$all writepsf protein_heme_mod.psf
$all delete
$cname delete
$hname delete
mol delete all

#mol new protein_heme_mod.psf
mol new protein_heme_mod.pdb
set protA [atomselect top "protein and segname PA"]
set protB [atomselect top "protein and segname PB"]
set hemA [atomselect top "resname HER"]
$protA writepdb protA.pdb
$protB writepdb protB.pdb
$hemA set resname HER
$hemA writepdb hemA.pdb

segment PA {
	pdb protA.pdb
	first ACE
	last CTER
}

segment HEMA {
	pdb hemA.pdb
	first none
	last none
}

segment PB {
        pdb protB.pdb
        first ACE
        last CTER
}
patch LCAB PA:67 HEMA:901
patch LCAC PA:70 HEMA:901
patch PHEM PA:71 HEMA:901
patch PHEM PA:110 HEMA:901

patch LCAB PA:100 HEMA:902
patch LCAC PA:103 HEMA:902
patch PHEM PA:104 HEMA:902
patch PHEM PA:167 HEMA:902

patch LCAB PA:136 HEMA:903
patch LCAC PA:139 HEMA:903
patch PHEM PA:140 HEMA:903
patch PHEM PA:85 HEMA:903

patch LCAB PA:160 HEMA:904
patch LCAC PA:163 HEMA:904
patch PHEM PA:164 HEMA:904
patch PHEM PA:216 HEMA:904

patch LCAB PA:183 HEMA:905
patch LCAC PA:186 HEMA:905
patch PHEM PA:187 HEMA:905
patch PHEM PA:153 HEMA:905

patch LCAB PA:209 HEMA:906
patch LCAC PA:212 HEMA:906
patch PHEM PA:213 HEMA:906
patch PHEM PA:262 HEMA:906

patch LCAB PA:233 HEMA:907
patch LCAC PA:236 HEMA:907
patch PHEM PA:237 HEMA:907
patch PHEM PA:200 HEMA:907

patch LCAB PA:255 HEMA:908
patch LCAC PA:258 HEMA:908
patch PHEM PA:259 HEMA:908
patch PHEM PA:321 HEMA:908

patch LCAB PA:278 HEMA:909
patch LCAC PA:281 HEMA:909
patch PHEM PA:282 HEMA:909
patch PHEM PA:248 HEMA:909

patch LCAB PA:313 HEMA:910
patch LCAC PA:316 HEMA:910
patch PHEM PA:287 HEMA:910
patch PHEM PA:317 HEMA:910

coordpdb protA.pdb PA
coordpdb hemA.pdb HEMA
coordpdb protB.pdb PB
guesscoord
writepdb protein_heme_mod2.pdb
writepsf protein_heme_mod2.psf
mol delete all
readpsf mem-wat-ion-ca.psf pdb mem-wat-ion-ca.pdb
readpsf protein_heme_mod2.psf pdb protein_heme_mod2.pdb
writepsf mem-prot_SO_reduced.psf
writepdb mem-prot_SO_reduced.pdb
quit
