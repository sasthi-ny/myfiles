package require psfgen
topology 	      ../../../toppar/top_all36_prot.rtf
topology 	      ../../../toppar/toppar_water_ions.str
topology 	      ../../../toppar/top_hemeOx_charmm36.inp
topology	      ../../../toppar/top_all36_lipid.rtf
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_moreions.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_ions_won.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_prot_na_combined.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_prot_retinol.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_sphingo.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_archaeal.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_bacterial.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_cardiolipin.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_dag.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_inositol.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_lnp.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_lps.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_mycobacterial.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_miscellaneous.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar/toppar_all36_lipid_model.str
topology              ../../../toppar_mem-pro-wat_charmmgui/toppar/toppar_all36_lipid_prot.str

# Aligning MtrCAB #
mol new protABC.psf
mol addfile protABC.pdb
set all1 [atomselect top all]
set protB1 [atomselect top "segname PB1 and name CA and resid 47 to 328 331 to 697"]
#$protB1 writepdb test1.pdb

# Load OM-MtrCAB(S. baltica) built from CHARMM-GUI
mol new OM-MtrCAB_S_baltica.psf
mol new OM-MtrCAB_S_baltica.pdb
set all2 [atomselect top all]
set protB2 [atomselect top "segname PROB and name CA"]
#$protB2 writepdb test2.pdb
set M1 [measure fit $protB1 $protB2]
$all1 move $M1
$all1 writepdb protABC_aligned.pdb
# End aligning #

mol new protABC_aligned.pdb
set chainA [atomselect top "segname PA"]
set chainB [atomselect top "segname PB1"]
set chainC [atomselect top "segname PC"]
set hem [atomselect top "resname HEC"]
set ca [atomselect top "resname CAL"]
set wat1 [atomselect top water]
$chainA writepdb protA.pdb
$chainB writepdb protB.pdb
$chainC writepdb protC.pdb
$hem writepdb heme.pdb
$ca writepdb ca.pdb
$wat1 writepdb wat.pdb
mol delete all

mol new OM-MtrCAB_S_baltica.psf
mol addfile OM-MtrCAB_S_baltica.pdb
#set all_not_prot [atomselect top "all not protein"]
set lpsa [atomselect top "chain L"]
set memb [atomselect top "segname MEMB"]
set wat2 [atomselect top water]
#set mg [atomselect top "resname MG"]
set ions [atomselect top "segname IONS"]
#set na [atomselect top "resname SOD"]
$all_not_prot writepdb step5_input_all_not_prot.pdb
$lpsa writepdb S_baltica_lpsa.pdb
$memb writepdb S_baltica_mem.pdb
$wat2 writepdb S_baltica_wat.pdb
$ions writepdb S_baltica_ions.pdb
mol delete all

#mol new protA.pdb
segment PA {pdb protA.pdb}
#mol new protB.pdb
segment PB {pdb protB.pdb}
#mol new protC.pdb
segment PC {pdb protC.pdb}
segment HEC {pdb heme.pdb}
segment CA {pdb ca.pdb}
segment WAT1 {pdb wat.pdb}
#mol new step5_lpsa.pdb
#segment LPSA {pdb step5_lpsa.pdb}
#mol new step5_mem.pdb
segment MEM {pdb S_baltica_mem.pdb}
#mol new step5_wat.pdb
segment WAT2 {pdb S_baltica_wat.pdb}
#mol new step5_mg.pdb
segment IONS {pdb S_baltica_ions.pdb}
#mol new step5_ca.pdb

coordpdb protA.pdb PA
coordpdb protB.pdb PB
coordpdb protC.pdb PC
coordpdb heme.pdb HEC
coordpdb ca.pdb CA
coordpdb wat.pdb WAT1
coordpdb S_baltica_lpsa.pdb LPSA
coordpdb S_baltica_mem.pdb MEM
coordpdb S_baltica_wat.pdb WAT2
coordpdb S_baltica_ions.pdb IONS
guesscoord
writepdb final_OM_MtrCAB.pdb
writepsf final_OM_MtrCAB.psf

quit
