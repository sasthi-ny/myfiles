package require psfgen

resetpsf
## Loading topology files
topology ../../../toppar/top_all36_prot.rtf
topology ../../../toppar/top_all36_prot_mod_res.rtf
topology ../../../toppar/toppar_water_ions.str
topology ../../../toppar/top_hemeOx_charmm36.inp

## Loading PDB file of MtrCAB of S. baltica
mol new 6r2q.pdb

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias residue CA CAL
pdbalias atom CA CA CAL
pdbalias residue HOH TIP3
pdbalias atom HOH O OH2

## Changing the resname of the residues covalently connected to heme groups
set protA1 [atomselect top "protein and chain A and resid 62 to 73"]
set cname [atomselect top "protein and chain A and resid 67 70"]
set hname [atomselect top "protein and chain A and resid 71"]
$cname set resname CYO
$hname set resname HSO
$protA1 writepdb protA_1.pdb
$protA1 delete
$cname delete
$hname delete

set protA2 [atomselect top "protein and chain A and resid 78 to 112"]
set cname [atomselect top "protein and chain A and resid 100 103"]
set hname [atomselect top "protein and chain A and resid 85 104 110"]
$cname set resname CYO
$hname set resname HSO
$protA2 writepdb protA_2.pdb
$protA2 delete
$cname delete
$hname delete

set protA3 [atomselect top "protein and chain A and resid 116 to 333"]
set cname [atomselect top "protein and chain A and resid 136 139 160 163 183 186 209 212 233 236 255 258 278 281 313 316"]
set hname [atomselect top "protein and chain A and resid 140 153 164 167 187 200 213 216 237 248 259 262 282 287 317 321"]
$cname set resname CYO
$hname set resname HSO
$protA3 writepdb protA_3.pdb
$protA3 delete

set watA [atomselect top "chain A and water"]
$watA writepdb watA.pdb
$watA delete

## Saving heme groups of chain A of SB
mol new 6r2q.pdb
set hem1 [atomselect top "resname HEC and chain A"]
#$hem1 set segname HEMA
$hem1 writepdb hemA.pdb
$hem1 delete
## End saving ##

# Aligning MtrA of alphafold with S.B. # The missing segment are taken from alphafold structure
mol new MtrA_from_alphafold.pdb
set all0 [atomselect top all]
set alpha1 [atomselect top "alpha and resid 62 to 73 78 to 112 116 to 333"]

mol new MtrA_prot_SB.pdb
set alpha2 [atomselect top "alpha and resid 62 to 73 78 to 112 116 to 333"]
set M1 [measure fit $alpha1 $alpha2]
$all0 move $M1
$all0 writepdb MtrA_from_alphafold_aligned.pdb
$all0 delete
$alpha1 delete
$alpha2 delete

# Saving new residues from alphafold structure which are to be added to MtrA of S.O. #
mol new MtrA_from_alphafold_aligned.pdb
set add_res [atomselect top "protein and chain A and resid 36 to 61"]
set add_res2 [atomselect top "protein and chain A and resid 74 to 77"]
set add_res3 [atomselect top "protein and chain A and resid 113 to 115"]
$add_res writepdb mtra_af_36_61.pdb
$add_res2 writepdb mtra_af_74_77.pdb
$add_res3 writepdb mtr_af_113_115.pdb

segment PA {
	pdb mtra_af_36_61.pdb
	pdb protA_1.pdb
	residue	74 SER A
	residue	75 GLU A
	residue	76 LYS A
	residue	77 VAL A
	pdb protA_2.pdb
	residue	113 GLY A
	residue	114 GLY A
	residue	115 ASN A
	pdb protA_3.pdb
	mutate 130 ASP
	mutate 146 MET
	mutate 151 GLY
	mutate 229 VAL
	mutate 231 ASP
	mutate 268 GLY
	first ACE
	last CTER
}

segment HEMA {
	pdb hemA.pdb
	first none
	last none
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
#regenerate angles dihedrals # No need

coordpdb mtra_af_36_61.pdb PA
coordpdb protA_1.pdb PA
coordpdb protA_2.pdb PA
coordpdb protA_3.pdb PA
coordpdb hemA.pdb HEMA
guesscoord

writepdb protAtemp.pdb
writepsf protAtemp.psf
mol delete all

# Aligning MtrB alphafold - MtrB S.B.#
mol new alpha_fold_mtrb.pdb
set all1 [atomselect top all]
set alpha1 [atomselect top "alpha and resid 47 to 328 331 to 697"]
$alpha1 set chain B
#$alpha1 writepdb test1.pdb
mol new 6r2q.pdb
set alpha2 [atomselect top "alpha and chain B and resid 47 to 328 329 to 695"]
#$alpha2 writepdb test2.pdb
set M1 [measure fit $alpha1 $alpha2]
$all1 move $M1
$all1 writepdb alpha_fold_mtrb_aligned.pdb
$all1 delete
$alpha1 delete
$alpha2 delete

mol delete all

mol new alpha_fold_mtrb_aligned.pdb
set resid_329_330 [atomselect top "resid 329 330"]
$resid_329_330 set chain B
$resid_329_330 writepdb MtrB_329_330.pdb
mol delete top

mol new 6r2q.pdb
set protB1 [atomselect top "protein and chain B and resid 47 to 328"]
$protB1 writepdb protB_1.pdb
$protB1 delete

set prot1 [atomselect top "protein and chain B and resid 329 to 695"]
$prot1 writepdb prot2B.pdb

set ca_ionsB [atomselect top "chain B and resname CA"]
#$ca_ions set "ATOM" "HETATM"
$ca_ionsB writepdb CA_ions.pdb
$ca_ionsB delete

set watB [atomselect top "chain B and (water and same residue as within 5 of resname CA)"]
$watB writepdb watB.pdb
$watB delete


# Open .pdb
mol new prot2B.pdb
set prot2 [atomselect top all]

# Get list of resids
set old_resid [ [ atomselect top all ] get resid ]

# Create a new list  {}
set new_resid {}
foreach r $old_resid {
  incr r 2
  lappend new_resid $r
}

# Reassign the resid using the new list
[ atomselect top all ] set resid $new_resid
# Write the PDB
[ atomselect top all ] writepdb prot2B_mod_resid.pdb

segment PB {
	pdb protB_1.pdb
	mutate	83	PHE
	mutate	87	LEU
	mutate	88	ASN
	mutate	90	LYS
	mutate	100	ASP
	mutate	112	ASP
	mutate	134	ASP
	mutate	140	SER
	mutate	143	ALA
	mutate	145	ILE
	mutate	147	GLY
	mutate	148	ASN
	mutate	157	ILE
	mutate	163	ASN
	mutate	174	ALA
	mutate	189	GLU
	mutate	199	TYR
	mutate	201	ASN
	mutate	212	GLN
	mutate	216	SER
	mutate	234	THR
	mutate	239	VAL
	mutate	245	ARG
	mutate	251	SER
	mutate	265	ASP
	mutate	267	GLU
	mutate	268	ASN
	mutate	279	GLN
	mutate	281	THR
	mutate	282	MET
	mutate	302	GLY
	mutate	303	SER
	mutate	305	ALA
	mutate	308	GLY
	mutate	310	ILE
	mutate	324	ASP
	mutate	325	ASN
	mutate	327	ARG
	pdb MtrB_329_330.pdb
	pdb prot2B_mod_resid.pdb
	mutate 332 LEU
	mutate 333 ASN
	mutate 335 ASP
	mutate 337 VAL
	mutate 344 LEU
	mutate 346 MET
	mutate 355 SER
	mutate 357 ASP
	mutate 361 THR
	mutate 367 TYR
	mutate 375 VAL
	mutate 381 ILE
	mutate 399 ARG
	mutate 400 THR
	mutate 414 ASP
	mutate 415 ILE
	mutate 424 LYS
	mutate 427 GLN
	mutate 429 ASP
	mutate 445 LEU
	mutate 451 ASP
	mutate 463 ASN
	mutate 468 GLN
	mutate 490 ASP
	mutate 500 ILE
	mutate 505 LEU
	mutate 512 VAL
	mutate 546 ALA
	mutate 551 THR
	mutate 558 THR
	mutate 581 GLN
	mutate 601 LEU
	mutate 609 ASN
	mutate 611 ASP
	mutate 646 LEU
	mutate 669 ASP
	mutate 679 SER
	mutate 689 LEU
	mutate 697 LEU
	first ACE
	last CTER
}
coordpdb protB_1.pdb PB
coordpdb MtrB_329_330.pdb PB
coordpdb prot2B_mod_resid.pdb PB
guesscoord
writepdb protBtemp.pdb

mol delete all

mol new CA_ions.pdb
segment CAB {
	pdb CA_ions.pdb
	first none
	last none
}
coordpdb CA_ions.pdb CAB
guesscoord

# Saving chain C of two mtrC #
mol new 6r2q.pdb
set mtrc_sb [atomselect top "protein and chain C"]
$mtrc_sb writepdb mtrc_sb_C.pdb
$mtrc_sb delete
mol new 4lm8_MtrC_SO.pdb
set mtrc_so [atomselect top all]
$mtrc_so set chain C
$mtrc_so writepdb mtrc_so_C.pdb
$mtrc_so delete

#set ca_so [atomselect top "resname CA"]
#$ca_so set chain C
#$ca_so set segname CAC
#$ca_so writepdb ca_so_c.pdb

mol delete all
# End saving #

# Aligning two MtrC #
mol new mtrc_so_C.pdb
set all1 [atomselect top all]
set alpha1 [atomselect top "alpha and resid 45 to 86 88 to 111 112 to 128 129 to 170 171 to 198 199 to 213 215 to 220 224 to 233 241 to 250 257 to 338 340 to 390 391 to 412 414 to 420 422 to 428 431 to 459 467 to 509 511 to 549 550 to 567 570 to 573 575 to 598 604 to 622 628 to 670"]

mol new mtrc_sb_C.pdb
set all2 [atomselect top all]
set alpha2 [atomselect top "alpha and resid 45 to 86 87 to 110 120 to 136 140 to 181 183 to 210 213 to 227 228 to 233 234 to 243 244 to 253 254 to 335 336 to 386 394 to 415 416 to 422 423 to 429 430 to 458 459 to 501 502 to 540 542 to 559 560 to 563 564 to 587 588 to 606 607 to 649"]

set M1 [measure fit $alpha1 $alpha2]
$all1 move $M1
#set crystal2 [$all1 move $M1] - not working
$all1 writepdb MtrC_SO_aligned.pdb
mol delete all
# Done aligning #

# Saving water #
mol new MtrC_SO_aligned.pdb
set watC [atomselect top "water and same residue as within 5 of resname CA"]
$watC writepdb watC.pdb
$watC delete
# End saving #


## Aligning heme groups ##
mol new 4lm8_MtrC_SO.pdb
set heme_so_all [atomselect top "resname HEC"]
set heme_so_fe [atomselect top "name FE and resname HEC"]
$heme_so_fe writepdb heme_so.pdb
mol new 6r2q.pdb
set heme_sb_fe [atomselect top "chain C and name FE and resname HEC"]
$heme_sb_fe writepdb heme_sb.pdb
set M2 [measure fit $heme_so_fe $heme_sb_fe]
$heme_so_all move $M2
$heme_so_all writepdb heme_so_c_aligned.pdb
$heme_so_fe delete
$heme_sb_fe delete
$heme_so_all delete
## Done aligning heme groups ##

mol delete all
mol new MtrC_SO_aligned.pdb
set prot [atomselect top "protein"]
set ca_ionsC [atomselect top "resname CA"] 
$prot writepdb MtrC_SO_aligned_prot.pdb
$ca_ionsC writepdb CA_ionsC.pdb
$prot delete
$ca_ionsC delete
mol delete top

mol new MtrC_SO_aligned_prot.pdb
set all [atomselect top all]
set cname [atomselect top "protein and resid 186 189 206 209 257 260 282 285 306 309 496 499 538 541 576 579 614 617 655 658"]
set hname [atomselect top "protein and resid 190 198 210 230 233 261 286 297 310 319 500 507 542 561 564 580 618 627 659 666"]
$cname set resname CYO
$hname set resname HSO
$all writepdb MtrC_SO_aligned_ch_res.pdb
segment PC {
	pdb MtrC_SO_aligned_ch_res.pdb
	first ACE
	last CT3
}

mol new heme_so_c_aligned.pdb
segment HEMC {
	pdb heme_so_c_aligned.pdb
	first none
	last none
}

## Adding patches between heme-prot ##
# Disulfide bond #
patch DISU PC:444 PC:453
# Connections with the heme groups #
patch LCAB PC:186 HEMC:801
patch LCAC PC:189 HEMC:801
patch PHEM PC:233 HEMC:801
patch PHEM PC:190 HEMC:801

#Added after visually seen for CAC
patch LCAB PC:206 HEMC:802
patch LCAC PC:209 HEMC:802  
patch PHEM PC:210 HEMC:802
patch PHEM PC:198 HEMC:802

patch LCAB PC:257 HEMC:803
patch LCAC PC:260 HEMC:803
patch PHEM PC:230 HEMC:803
patch PHEM PC:261 HEMC:803

patch LCAB PC:282 HEMC:804
patch LCAC PC:285 HEMC:804
patch PHEM PC:286 HEMC:804
patch PHEM PC:319 HEMC:804

#Added from visually seen for CAC
patch LCAB PC:306 HEMC:805
patch LCAC PC:309 HEMC:805 
patch PHEM PC:310 HEMC:805
patch PHEM PC:297 HEMC:805

#Added from visually seen for CAC
patch LCAB PC:496 HEMC:806
patch LCAC PC:499 HEMC:806 
patch PHEM PC:500 HEMC:806
patch PHEM PC:564 HEMC:806

patch LCAB PC:538 HEMC:807
patch LCAC PC:541 HEMC:807
patch PHEM PC:542 HEMC:807
patch PHEM PC:507 HEMC:807

patch LCAB PC:576 HEMC:808
patch LCAC PC:579 HEMC:808
patch PHEM PC:561 HEMC:808
patch PHEM PC:580 HEMC:808

patch LCAB PC:614 HEMC:809
patch LCAC PC:617 HEMC:809
patch PHEM PC:618 HEMC:809
patch PHEM PC:666 HEMC:809

patch LCAB PC:655 HEMC:810
patch LCAC PC:658 HEMC:810
patch PHEM PC:627 HEMC:810
patch PHEM PC:659 HEMC:810
#regenerate angles dihedrals
# End adding patches #
coordpdb MtrC_SO_aligned_ch_res.pdb PC
coordpdb heme_so_c_aligned.pdb HEMC
guesscoord

mol new CA_ionsC.pdb
segment CAC {
	pdb CA_ionsC.pdb
	first none
	last none
}
#mol new watA.pdb - no need
segment WA {
	pdb watA.pdb
	first none
	last none
}
segment WB {
	pdb watB.pdb
	first none
	last none
}
segment WC {
	pdb watC.pdb
	first none
	last none
}
coordpdb CA_ionsC.pdb CAC
coordpdb watA.pdb WA
coordpdb watB.pdb WB
coordpdb watC.pdb WC
guesscoord
writepdb protABC.pdb
writepsf protABC.psf

## Removing water from protABC.pdb ##
mol new protABC.psf
mol addfile protABC.pdb
#set all [atomselect top all]
#$all set beta 0
#$all writepdb protABC_beta0.pdb
set all_nowat [atomselect top "all not water"]
$all_nowat writepdb protABC_nowat.pdb
$all_nowat writepsf protABC_nowat.psf
$all_nowat delete
$all delete
# set beta 1 for the residues which are restrained #
mol new protABC.psf
mol addfile protABC.pdb
set all [atomselect top all]
$all set beta 1
$all writepdb protABC_beta1.pdb
$all delete

mol new protABC.psf
mol addfile protABC_beta1.pdb
set all [atomselect top all]
set mod1 [atomselect top "segname PA and resid 36 to 62 73 to 78 112 to 116 130 146 151 229 231 268"]
set mod2 [atomselect top "segname PB and resid 83 87 88 90 100 112 134 140 143 145 147 148 157 163 174 189 199 201 212 216 234 239 245 251 265 267 268 279 281  282 302 303 305 308 310 324 325 327 328 329 330 331 332 333 335 337 344 346 355 357 361 367 375 381 399 400 414 415 424 427 429 445 451 463 468 490 500 505 512 546 551 558 581 601 609 611 646 669 679 689 697"]
$mod1 set beta 0
$mod2 set beta 0
$all writepdb protABC_beta0_1.pdb
$all writepsf protABC_beta0_1.psf
$all delete
mol delete all

#mol new protABC_beta0_1.pdb
#coordpdb protABC_beta0_1.pdb
#guesscoord
#writepsf protABC_beta0_1.psf
#mol delete all

file delete CA_ions.pdb
file delete protA_1.pdb protA_2.pdb protA_3.pdb CA_ionsC.pdb protAtemp.pdb watA.pdb hemA.pdb
file delete prot2B_mod_resid.pdb mtra_af_74_77.pdb mtr_af_113_115.pdb protB2temp.pdb
file delete protB_1.pdb protBtemp.pdb watB.pdb MtrB_329_330.pdb prot2B.pdb CA_ionsB.pdb ca_so_c.pdb
file delete heme_so.pdb heme_sb.pdb heme_so_c_aligned.pdb
## End removing water ##
quit
