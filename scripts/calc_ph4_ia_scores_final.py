# -*- mode: python; tab-width: 4 -*- 

## 
# calc_ph4_ia_scores_red.py 
#
# Copyright (C) 2022 Thomas A. Seidel <thomas.seidel@univie.ac.at>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; see the file COPYING. If not, write to
# the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
##

import argparse
import os
import sys
import math

import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Pharm as Pharm
import CDPL.Math as Math
import CDPL.MolProp as MolProp
import CDPL.ForceField as ForceField
import CDPL.ConfGen as ConfGen


LIG_ENV_MAX_RADIUS = 20.0
REMOVE_NON_STD_RESIDUES = True
NON_STD_RESIDUE_MAX_ATOM_COUNT = 8
FTR_DIST_CUTOFF = 10.0
VDW_DIST_CUTOFF = 20.0
ESTAT_DIST_CUTOFF = 20.0
DIELECTRIC_CONST = 1.0
HYDROPHOBICIY_THRESH = 0.15


# AT -> (eps, Rmin/2)
CHARMM_LJ_PARAMS = {
    'C' :   (-0.110000, 2.000000),
    'CA' :  (-0.070000, 1.992400),
    'CC' :  (-0.070000, 2.000000),
    'CD' :  (-0.070000, 2.000000),
    'CE1' : (-0.068000, 2.090000),  
    'CE2' : (-0.064000, 2.080000),  
    'CP1' : (-0.020000, 2.275000),
    'CP2' : (-0.055000, 2.175000),
    'CP3' : (-0.055000, 2.175000),
    'CPH1' :(-0.050000, 1.800000),
    'CPH2' :(-0.050000, 1.800000),
    'CS' :  (-0.110000, 2.200000),
    'CPT' : (-0.099000, 1.860000),
    'CY' :  (-0.073000, 1.990000),
    'CAI' : (-0.073000, 1.990000),
    'CT' :  (-0.0200,   2.275),  
    'CT1' : (-0.0320,   2.000),
    'CT2' : (-0.0560,   2.010),
    'CT2A' :(-0.0560,   2.010),
    'CT3' : (-0.0780,   2.040),
    'H' :   (-0.046000, 0.224500),
    'HA' :  (-0.022000, 1.320000),
    'HB1' : (-0.022000, 1.320000),  
    'HB2' : (-0.028000, 1.340000),  
    'HE1' : (-0.031000, 1.250000),  
    'HE2' : (-0.026000, 1.260000),  
    'HC' :  (-0.046000, 0.224500),
    'HP' :  (-0.030000, 1.358200),
    'HR1' : (-0.046000, 0.900000),
    'HR2' : (-0.046000, 0.700000),
    'HR3' : (-0.007800, 1.468000),
    'HS' :  (-0.100000, 0.450000),
    'HA1' : (-0.045,    1.3400),
    'HA2' : (-0.034,    1.3400),
    'HA3' : (-0.024,    1.3400),
    'N' :   (-0.200000, 1.850000),
    'NC2' : (-0.200000, 1.850000),
    'NH1' : (-0.200000, 1.850000),
    'NH2' : (-0.200000, 1.850000),
    'NH3' : (-0.200000, 1.850000),
    'NP' :  (-0.200000, 1.850000),
    'NR1' : (-0.200000, 1.850000),
    'NR2' : (-0.200000, 1.850000),
    'NR3' : (-0.200000, 1.850000),
    'NY' :  (-0.200000, 1.850000),
    'O' :   (-0.120000, 1.700000),
    'OB' :  (-0.120000, 1.700000),
    'OC' :  (-0.120000, 1.700000),
    'OH1' : (-0.152100, 1.770000),
    'OS' :  (-0.152100, 1.770000),
    'S' :   (-0.450000, 2.000000),
    'SM' :  (-0.380000, 1.975000),
    'SS' :  (-0.470000, 2.200000) }

Pharm.FeatureType.HBD_N3 = Pharm.FeatureType.MAX_TYPE + 1
Pharm.FeatureType.HBD_N2 = Pharm.FeatureType.MAX_TYPE + 2
Pharm.FeatureType.HBD_Nar = Pharm.FeatureType.MAX_TYPE + 3
Pharm.FeatureType.HBD_Nam = Pharm.FeatureType.MAX_TYPE + 4
Pharm.FeatureType.HBD_Npl3 = Pharm.FeatureType.MAX_TYPE + 5
Pharm.FeatureType.HBD_N4 = Pharm.FeatureType.MAX_TYPE + 6
Pharm.FeatureType.HBD_O3 = Pharm.FeatureType.MAX_TYPE + 7
Pharm.FeatureType.HBD_S3 = Pharm.FeatureType.MAX_TYPE + 8
Pharm.FeatureType.HBA_N3 = Pharm.FeatureType.MAX_TYPE + 9
Pharm.FeatureType.HBA_N2 = Pharm.FeatureType.MAX_TYPE + 10
Pharm.FeatureType.HBA_N1 = Pharm.FeatureType.MAX_TYPE + 11
Pharm.FeatureType.HBA_Nar = Pharm.FeatureType.MAX_TYPE + 12
Pharm.FeatureType.HBA_Npl3 = Pharm.FeatureType.MAX_TYPE + 13
Pharm.FeatureType.HBA_O3 = Pharm.FeatureType.MAX_TYPE + 14
Pharm.FeatureType.HBA_O2 = Pharm.FeatureType.MAX_TYPE + 15
Pharm.FeatureType.HBA_Oco2 = Pharm.FeatureType.MAX_TYPE + 16
Pharm.FeatureType.HBA_S3 = Pharm.FeatureType.MAX_TYPE + 17
Pharm.FeatureType.HBA_S2 = Pharm.FeatureType.MAX_TYPE + 18
Pharm.FeatureType.HBD_N = Pharm.FeatureType.MAX_TYPE + 19
Pharm.FeatureType.HBD_O = Pharm.FeatureType.MAX_TYPE + 20
Pharm.FeatureType.HBD_S = Pharm.FeatureType.MAX_TYPE + 21
Pharm.FeatureType.HBA_N = Pharm.FeatureType.MAX_TYPE + 22
Pharm.FeatureType.HBA_O = Pharm.FeatureType.MAX_TYPE + 23
Pharm.FeatureType.HBA_S = Pharm.FeatureType.MAX_TYPE + 24


FEATURE_NAMES = { Pharm.FeatureType.POSITIVE_IONIZABLE : 'PI',
                  Pharm.FeatureType.NEGATIVE_IONIZABLE : 'NI',
                  Pharm.FeatureType.AROMATIC : 'AR',
                  Pharm.FeatureType.HYDROPHOBIC : 'H',
                  Pharm.FeatureType.H_BOND_DONOR : 'HBD',
                  Pharm.FeatureType.H_BOND_ACCEPTOR : 'HBA',
                  Pharm.FeatureType.HALOGEN_BOND_DONOR : 'XBD',
                  Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR : 'XBA',
                  Pharm.FeatureType.HBD_N : 'HBD_N',
                  Pharm.FeatureType.HBD_O : 'HBD_O',
                  Pharm.FeatureType.HBD_S : 'HBD_S',
                  Pharm.FeatureType.HBA_N : 'HBA_N',
                  Pharm.FeatureType.HBA_O : 'HBA_O',
                  Pharm.FeatureType.HBA_S : 'HBA_S',
                  Pharm.FeatureType.HBD_N3 : 'HBD_N3',
                  Pharm.FeatureType.HBD_N2 : 'HBD_N2',
                  Pharm.FeatureType.HBD_Nar : 'HBD_Nar',
                  Pharm.FeatureType.HBD_Nam : 'HBD_am',
                  Pharm.FeatureType.HBD_Npl3 : 'HBD_Npl3',
                  Pharm.FeatureType.HBD_N4 : 'HBD_N4',
                  Pharm.FeatureType.HBD_O3 : 'HBD_O3',
                  Pharm.FeatureType.HBD_S3 : 'HBD_S3',
                  Pharm.FeatureType.HBA_N3 : 'HBA_N3',
                  Pharm.FeatureType.HBA_N2 : 'HBA_N2',
                  Pharm.FeatureType.HBA_N1 : 'HBA_N1',
                  Pharm.FeatureType.HBA_Nar : 'HBA_Nar',
                  Pharm.FeatureType.HBA_Npl3 : 'HBA_Npl3',
                  Pharm.FeatureType.HBA_O3 : 'HBA_O3',
                  Pharm.FeatureType.HBA_O2 : 'HBA_O2',
                  Pharm.FeatureType.HBA_Oco2 : 'HBA_Oco2',
                  Pharm.FeatureType.HBA_S3 : 'HBA_S3',
                  Pharm.FeatureType.HBA_S2 : 'HBA_S2' }

LIG_HBA_FEATURE_TYPES = [ Pharm.FeatureType.H_BOND_ACCEPTOR,
                          Pharm.FeatureType.HBA_N3,
                          Pharm.FeatureType.HBA_N2,
                          Pharm.FeatureType.HBA_N1,
                          Pharm.FeatureType.HBA_Nar,
                          Pharm.FeatureType.HBA_Npl3,
                          Pharm.FeatureType.HBA_O3,
                          Pharm.FeatureType.HBA_O2,
                          Pharm.FeatureType.HBA_Oco2,
                          Pharm.FeatureType.HBA_S3,
                          Pharm.FeatureType.HBA_S2 ]

LIG_HBD_FEATURE_TYPES = [ Pharm.FeatureType.H_BOND_DONOR,
                          Pharm.FeatureType.HBD_N3,
                          Pharm.FeatureType.HBD_N2,
                          Pharm.FeatureType.HBD_Nar,
                          Pharm.FeatureType.HBD_Nam,
                          Pharm.FeatureType.HBD_Npl3,
                          Pharm.FeatureType.HBD_N4,
                          Pharm.FeatureType.HBD_O3,
                          Pharm.FeatureType.HBD_S3 ]

ENV_HBA_FEATURE_TYPES = [ Pharm.FeatureType.HBA_N,
                          Pharm.FeatureType.HBA_O,
                          Pharm.FeatureType.HBA_S ]

ENV_HBD_FEATURE_TYPES = [ Pharm.FeatureType.HBD_N,
                          Pharm.FeatureType.HBD_O,
                          Pharm.FeatureType.HBD_S ]
                          
SYBYL_TO_HBD_TYPE = { Chem.SybylAtomType.N_3 : Pharm.FeatureType.HBD_N3,
                      Chem.SybylAtomType.N_2 : Pharm.FeatureType.HBD_N2,
                      Chem.SybylAtomType.N_ar : Pharm.FeatureType.HBD_Nar,
                      Chem.SybylAtomType.N_am : Pharm.FeatureType.HBD_Nam,
                      Chem.SybylAtomType.N_pl3 : Pharm.FeatureType.HBD_Npl3,
                      Chem.SybylAtomType.N_4 : Pharm.FeatureType.HBD_N4,
                      Chem.SybylAtomType.O_3 : Pharm.FeatureType.HBD_O3,
                      Chem.SybylAtomType.S_3 : Pharm.FeatureType.HBD_S3 }

SYBYL_TO_HBA_TYPE = { Chem.SybylAtomType.N_3 : Pharm.FeatureType.HBA_N3,
                      Chem.SybylAtomType.N_2 : Pharm.FeatureType.HBA_N2,
                      Chem.SybylAtomType.N_1 : Pharm.FeatureType.HBA_N1,
                      Chem.SybylAtomType.N_ar : Pharm.FeatureType.HBA_Nar,
                      Chem.SybylAtomType.N_pl3 : Pharm.FeatureType.HBA_Npl3,
                      Chem.SybylAtomType.O_3 : Pharm.FeatureType.HBA_O3,
                      Chem.SybylAtomType.O_2 : Pharm.FeatureType.HBA_O2,
                      Chem.SybylAtomType.O_co2 : Pharm.FeatureType.HBA_Oco2,
                      Chem.SybylAtomType.S_3 : Pharm.FeatureType.HBA_S3,
                      Chem.SybylAtomType.S_2 : Pharm.FeatureType.HBA_S2 }

ATOM_TO_HBD_TYPE = { Chem.AtomType.N : Pharm.FeatureType.HBD_N,
                     Chem.AtomType.O : Pharm.FeatureType.HBD_O,
                     Chem.AtomType.S : Pharm.FeatureType.HBD_S }

ATOM_TO_HBA_TYPE = { Chem.AtomType.N : Pharm.FeatureType.HBA_N,
                     Chem.AtomType.O : Pharm.FeatureType.HBA_O,
                     Chem.AtomType.S : Pharm.FeatureType.HBA_S }

INTERACTION_TYPES = [ (Pharm.FeatureType.POSITIVE_IONIZABLE, Pharm.FeatureType.AROMATIC, Pharm.CationPiInteractionScore(False), True),
                      (Pharm.FeatureType.AROMATIC, Pharm.FeatureType.POSITIVE_IONIZABLE, Pharm.CationPiInteractionScore(True), True),
                      (Pharm.FeatureType.HYDROPHOBIC, Pharm.FeatureType.HYDROPHOBIC, Pharm.HydrophobicInteractionScore(), False),
                      (Pharm.FeatureType.AROMATIC, Pharm.FeatureType.AROMATIC, Pharm.FeatureInteractionScoreCombiner(Pharm.OrthogonalPiPiInteractionScore(), Pharm.ParallelPiPiInteractionScore()), True),

                      (Pharm.FeatureType.H_BOND_DONOR, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.H_BOND_DONOR, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.H_BOND_DONOR, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),

                      (Pharm.FeatureType.HBD_N3, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_N3, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_N3, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),

                      (Pharm.FeatureType.HBD_N2, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_N2, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_N2, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),
                      
                      (Pharm.FeatureType.HBD_Nar, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_Nar, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_Nar, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),

                      (Pharm.FeatureType.HBD_Nam, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_Nam, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_Nam, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),

                      (Pharm.FeatureType.HBD_Npl3, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_Npl3, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_Npl3, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),

                      (Pharm.FeatureType.HBD_N4, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_N4, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_N4, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),
                
                      (Pharm.FeatureType.HBD_O3, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_O3, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_O3, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),

                      (Pharm.FeatureType.HBD_S3, Pharm.FeatureType.HBA_N, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_S3, Pharm.FeatureType.HBA_O, Pharm.HBondingInteractionScore(True), False),
                      (Pharm.FeatureType.HBD_S3, Pharm.FeatureType.HBA_S, Pharm.HBondingInteractionScore(True), False),

                      (Pharm.FeatureType.H_BOND_ACCEPTOR, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.H_BOND_ACCEPTOR, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.H_BOND_ACCEPTOR, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                      
                      (Pharm.FeatureType.HBA_N3, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_N3, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_N3, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                       
                      (Pharm.FeatureType.HBA_N2, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_N2, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_N2, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                        
                      (Pharm.FeatureType.HBA_N1, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_N1, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_N1, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                      
                      (Pharm.FeatureType.HBA_Nar, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_Nar, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_Nar, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                        
                      (Pharm.FeatureType.HBA_Npl3, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_Npl3, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_Npl3, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                      
                      (Pharm.FeatureType.HBA_O3, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_O3, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_O3, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                     
                      (Pharm.FeatureType.HBA_O2, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_O2, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_O2, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                    
                      (Pharm.FeatureType.HBA_Oco2, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_Oco2, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_Oco2, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                    
                      (Pharm.FeatureType.HBA_S3, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_S3, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_S3, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                       
                      (Pharm.FeatureType.HBA_S2, Pharm.FeatureType.HBD_N, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_S2, Pharm.FeatureType.HBD_O, Pharm.HBondingInteractionScore(False), False),
                      (Pharm.FeatureType.HBA_S2, Pharm.FeatureType.HBD_S, Pharm.HBondingInteractionScore(False), False),
                     
                      (Pharm.FeatureType.HALOGEN_BOND_DONOR, Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR, Pharm.XBondingInteractionScore(True), True) ]

def parseArguments():
    parser = argparse.ArgumentParser(description='Calculates GRAIL interactions scores between ligand and environment pharmacophore features for a set of input ligand-protein complexes.')
    parser.add_argument('-d',
                        dest='complex_data_dir',
                        required=True,
                        help='[Required] The directory containing the ligand-protein complexes to process, organized in PDBBind manner.',
                        nargs=1)
    parser.add_argument('-o',
                        dest='out_csv_file',
                        required=True,
                        help='[Required] The path of the output CSV-file containing the interaction scores calculated for each input complex.',
                        nargs=1)
    parser.add_argument('-p',
                        dest='output_ph4s',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Store generated pharmacophores.')
    parser.add_argument('-x',
                        dest='exact',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Regard ligand-side feature orientation when calculating feature interaction scores.')
    parser.add_argument('-r',
                        dest='dyn_hbd',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Generate dynamic HBD features on het-atom H-rotors.')

    parse_args = parser.parse_args()
    return parse_args

def calcScoresForInteraction(ia_type, lig_ph4, env_ph4, exact):
    score_sum = 0.0
    score_max_sum = 0.0

    for lig_ftr in lig_ph4:
        if Pharm.getType(lig_ftr) != ia_type[0]:
            continue

        max_score = 0.0
        lig_ftr_pos = Chem.get3DCoordinates(lig_ftr)
        lig_ftr_wt = Pharm.getWeight(lig_ftr)
        
        for env_ftr in env_ph4:
            if Pharm.getType(env_ftr) != ia_type[1]:
                continue

            env_ftr_pos = Chem.get3DCoordinates(env_ftr)

            if Math.length(env_ftr_pos - lig_ftr_pos) > FTR_DIST_CUTOFF:
                continue

            if exact:
                score = lig_ftr_wt * ia_type[2](lig_ftr, env_ftr)
            else:
                score = lig_ftr_wt * ia_type[2](lig_ftr_pos, env_ftr)
            
            score_sum += score
            max_score = max(score, max_score)

        score_max_sum += max_score
            
    return (score_sum, score_max_sum)

def calcEnvHDonorAcceptorOccupation(lig, ftr_type, env_ph4, acc_ftr):
    score_sum = 0.0
    score_max_sum = 0.0
    scoring_func = Pharm.HBondingInteractionScore(acc_ftr)

    for atom in lig.atoms:
        if Chem.getType(atom) == Chem.AtomType.H:
            continue

        atom_pos = Chem.get3DCoordinates(atom)
        max_score = 0.0
        
        for env_ftr in env_ph4:
            if Pharm.getType(env_ftr) != ftr_type:
                continue

            env_ftr_pos = Chem.get3DCoordinates(env_ftr)

            if Math.length(env_ftr_pos - atom_pos) > FTR_DIST_CUTOFF:
                continue

            score = scoring_func(atom_pos, env_ftr)
            
            score_sum += score
            max_score = max(score, max_score)

        score_max_sum += max_score
            
    return (score_sum, score_max_sum)

def calcElectrostaticInteractionEnergy(ligand, lig_env):
    ForceField.perceiveMMFF94AromaticRings(ligand, False)
    ForceField.perceiveMMFF94AromaticRings(lig_env, False)
    ForceField.assignMMFF94AtomTypes(ligand, False, False)
    ForceField.assignMMFF94AtomTypes(lig_env, False, False)
    ForceField.assignMMFF94BondTypeIndices(ligand, False, False)
    ForceField.assignMMFF94BondTypeIndices(lig_env, False, False)
    ForceField.calcMMFF94AtomCharges(ligand, False, False)
    ForceField.calcMMFF94AtomCharges(lig_env, False, False)

    energy = 0.0
    
    for lig_atom in ligand.atoms:
        lig_atom_pos = Chem.get3DCoordinates(lig_atom)
        lig_atom_charge = ForceField.getMMFF94Charge(lig_atom)

        if abs(lig_atom_charge) < 0.01:
            continue
        
        for env_atom in lig_env.atoms:
            env_atom_charge = ForceField.getMMFF94Charge(env_atom)

            if env_atom_charge == 0.0:
                continue
            
            env_atom_pos = Chem.get3DCoordinates(env_atom)
            dist = Math.length(env_atom_pos - lig_atom_pos)

            if dist > ESTAT_DIST_CUTOFF:
                continue

            energy += lig_atom_charge * env_atom_charge / (dist * DIELECTRIC_CONST)

    return energy

def getCHARMMAtomType(atom, molgraph, is_ligand):
    res_atom_name = ''
    res_code = ''

    if not is_ligand:
        if Biomol.hasResidueAtomName(atom):
            res_atom_name = Biomol.getResidueAtomName(atom)
        if Biomol.hasResidueCode(atom):
            res_code = Biomol.getResidueCode(atom)
    
    atomic_no = Chem.getType(atom)

    if atomic_no == Chem.AtomType.H:
        if not is_ligand:
            if res_code == 'HIS':
                if res_atom_name == 'HG' or res_atom_name == 'HD2':
                    return 'HR1'

                if res_atom_name == 'HE1':
                    return 'HR2'
                
        if atom.getNumAtoms() == 1 and MolProp.getBondCount(atom, molgraph) == 1:
            con_atom = atom.getAtom(0)

            if Chem.getAromaticityFlag(con_atom):
                return 'HP'

            con_atom_type = Chem.getType(con_atom)

            if con_atom_type == Chem.AtomType.S:
                return 'HS'

            if con_atom_type == Chem.AtomType.O:
                return 'H'

            if con_atom_type == Chem.AtomType.N:
                if not is_ligand and 'HT' in res_atom_name:
                    return 'HC'
                
                return 'H'
            
            if con_atom_type == Chem.AtomType.C:
                 h_count = MolProp.getAtomCount(con_atom, molgraph, Chem.AtomType.H)
                 hyb_state = Chem.getHybridizationState(con_atom)
                 
                 if hyb_state == Chem.HybridizationState.SP3:
                     if h_count == 3:
                         return 'HA3'

                     if h_count == 2:
                         return 'HA2'

                     if h_count == 1:
                         return 'HA1'
                     
                 elif hyb_state == Chem.HybridizationState.SP2:
                     if h_count == 2:
                         return 'HE2'

                     if h_count == 1:
                         return 'HE1'
                     
        return 'HA'

    if atomic_no == Chem.AtomType.C:
        if not is_ligand:
            if res_code == 'HIS':
                if res_atom_name == 'CG' or res_atom_name == 'CD2':
                    return 'CPH1'

                if res_atom_name == 'CE1':
                    return 'CPH2'
                
            elif res_code == 'PRO':
                if res_atom_name == 'CA':
                    return 'CP1'

                if res_atom_name == 'CB' or res_atom_name == 'CG':
                    return 'CP2'
                
                if res_atom_name == 'CD':
                    return 'CP3'

        if Chem.getAromaticityFlag(atom):
            if not is_ligand and res_code == 'TRP':
                if MolProp.getNumContainingSSSRRings(atom, molgraph) == 2:
                    return 'CPT'
            
            return 'CA'

        if MolProp.isCarbonylLikeAtom(atom, molgraph):
            if not is_ligand and res_code in [ 'ASN', 'ASP', 'GLN', 'GLU' ]:
                return 'CC'
                
            return 'C'
        
        if Chem.getHybridizationState(atom) == Chem.HybridizationState.SP3:
            if MolProp.getBondCount(atom, molgraph, 1, Chem.AtomType.S) == 1:
                return 'CS'
            
            h_count = MolProp.getAtomCount(atom, molgraph, Chem.AtomType.H)

            if h_count == 3:
                return 'CT3'

            if h_count == 2:
                return 'CT2'

            if h_count == 1:
                return 'CT1'

            if h_count == 0:
                return 'CT'

        elif Chem.getHybridizationState(atom) == Chem.HybridizationState.SP2:
            h_count = MolProp.getAtomCount(atom, molgraph, Chem.AtomType.H)

            if h_count == 1:
                return 'CE1'

            if h_count == 2:
                return 'CE2'
        
        return 'CT2'

    if atomic_no == Chem.AtomType.N:
        return 'N'

    if atomic_no == Chem.AtomType.O:
        if MolProp.getBondCount(atom, molgraph, 2, Chem.AtomType.C) == 1:
            return 'O'
            
        return 'OH1'

    if atomic_no == Chem.AtomType.S:
        if MolProp.getBondCount(atom, molgraph, 1, Chem.AtomType.C) == 1 and MolProp.getBondCount(atom, molgraph, 1, Chem.AtomType.S) == 1:
            return 'SM'

        if MolProp.getBondCount(atom, molgraph, 1) == 2 and MolProp.getBondCount(atom, molgraph, 1, Chem.AtomType.H) >= 1:
            return 'SS'
        
        return 'S'

    return ''
    
def getCHARMMLJParameters(molgraph, is_ligand):
    params = []

    for atom in molgraph.atoms:
        atom_type = getCHARMMAtomType(atom, molgraph, is_ligand)

        if atom_type == '':
            params.append((0.0, 0.0))
        else:
            params.append(CHARMM_LJ_PARAMS[atom_type])
    
    return params

# Uses CHARMM formula:
# V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2 * (Rmin,i,j/ri,j)**6]
# Eps,i,j = sqrt(eps,i * eps,j)
# Rmin,i,j = Rmin/2,i + Rmin/2,j
def calcLJVdWInteractionEnergy(ligand, lig_env):
    lig_atom_params = getCHARMMLJParameters(ligand, True)
    env_atom_params = getCHARMMLJParameters(lig_env, False)

    energy_att = 0.0
    energy_rep = 0.0

    for i in range(0, ligand.numAtoms):
        lig_atom_pos = Chem.get3DCoordinates(ligand.getAtom(i))
        lig_atom_param = lig_atom_params[i]

        for j in range(0, lig_env.numAtoms):
            env_atom_pos = Chem.get3DCoordinates(lig_env.getAtom(j))
            r_ij = Math.length(env_atom_pos - lig_atom_pos)

            if r_ij > VDW_DIST_CUTOFF:
                continue
    
            rel_dist_rep = (lig_atom_params[i][1] + env_atom_params[j][1]) / r_ij
            rel_dist_att = (lig_atom_params[i][1] + env_atom_params[j][1]) / r_ij
            
            eps = math.sqrt(lig_atom_params[i][0] * env_atom_params[j][0])
            energy_att += -eps * 2 * (rel_dist_att ** 6)
            energy_rep += eps * (rel_dist_rep ** 12)

    return (energy_att + energy_rep)

# Uses a Morse potential:
# V = Eps,i,j * (e^(-2 * alpha * (ri,j - Rmin,i,j)) - 2 * e^(-alpha * (ri,j - Rmin,i,j)))
# Eps,i,j = sqrt(eps,i * eps,j)
# Rmin,i,j = (Rmin/2,i + Rmin/2,j) * 2^(1/6)
def calcMorseVdWInteractionEnergy(ligand, lig_env):
    rmin_fact = math.pow(2, 1.0 / 6)
    alpha = 1.1

    lig_atom_params = getCHARMMLJParameters(ligand, True)
    env_atom_params = getCHARMMLJParameters(lig_env, False)
    energy = 0.0

    for i in range(0, ligand.numAtoms):
        lig_atom_pos = Chem.get3DCoordinates(ligand.getAtom(i))
        lig_atom_param = lig_atom_params[i]

        for j in range(0, lig_env.numAtoms):
            env_atom_pos = Chem.get3DCoordinates(lig_env.getAtom(j))
            r_ij = Math.length(env_atom_pos - lig_atom_pos)

            if r_ij > VDW_DIST_CUTOFF:
                continue
    
            r_min = rmin_fact * (lig_atom_params[i][1] + env_atom_params[j][1])
            r_delta = r_ij - r_min
            eps = math.sqrt(lig_atom_params[i][0] * env_atom_params[j][0])

            energy += eps * (math.exp(-2 * alpha * r_delta) - 2 * math.exp(-alpha * r_delta))

    return energy

def removeNonStdResidues(pdb_code, lig_env):
    residues = Biomol.ResidueList(lig_env)
    
    for res in residues:
        is_std_res = Biomol.ResidueDictionary.isStdResidue(Biomol.getResidueCode(res))
        
        if is_std_res and res.numAtoms < 5:
            print('!! While processing complex %s: isolated standard residue fragment of size %s found' % (pdb_code, str(res.numAtoms)), file=sys.stderr)
            lig_env -= res
            
        elif REMOVE_NON_STD_RESIDUES and not is_std_res and res.numAtoms <= NON_STD_RESIDUE_MAX_ATOM_COUNT:
            if res.numAtoms == 1 and Chem.AtomDictionary.isMetal(Chem.getType(res.atoms[0])):
                continue

            lig_env -= res

def checkProtein(pdb_code, protein):
    for atom in protein.atoms:
        if Chem.getType(atom) == Chem.AtomType.H and atom.numAtoms == 0:
            print('!! While processing complex %s: isolated hydrogen atom encountered' % pdb_code, file=sys.stderr)
             
        elif Chem.getType(atom) == Chem.AtomType.UNKNOWN:
            print('!! While processing complex %s: atom of unknown element encountered' % pdb_code, file=sys.stderr)

def prepareMolecule(mol, is_lig):
    Chem.ProtonationStateStandardizer().standardize(mol, Chem.ProtonationStateStandardizer.PHYSIOLOGICAL_CONDITION_STATE)
    Pharm.prepareForPharmacophoreGeneration(mol)

    if is_lig:
        Chem.perceiveSybylAtomTypes(mol, True)
    
    logp_calc = MolProp.XLogPCalculator()
    logp = logp_calc.calculate(mol)
    atom_hyds = logp_calc.getAtomContributions()

    for i in range(0, mol.numAtoms):
        MolProp.setHydrophobicity(mol.getAtom(i), atom_hyds(i))

    return logp

def generatePharmacophore(mol, ph4, dyn_hbd, is_lig):
    if not dyn_hbd:
        ph4_gen = Pharm.DefaultPharmacophoreGenerator(Pharm.DefaultPharmacophoreGenerator.Configuration(Pharm.DefaultPharmacophoreGenerator.Configuration.STATIC_H_DONORS | Pharm.DefaultPharmacophoreGenerator.Configuration.PI_NI_ON_CHARGED_GROUPS_ONLY))
    else:
        ph4_gen = Pharm.DefaultPharmacophoreGenerator(Pharm.DefaultPharmacophoreGenerator.Configuration.PI_NI_ON_CHARGED_GROUPS_ONLY)

    h_gen = Pharm.HydrophobicAtomFeatureGenerator()
    h_gen.setHydrophobicityThreshold(HYDROPHOBICIY_THRESH)

    ph4_gen.setFeatureGenerator(Pharm.FeatureType.HYDROPHOBIC, h_gen)

    if not is_lig:
        ph4_gen.enableFeature(Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR, True);

    ph4_gen.generate(mol, ph4)

def getFeatureAtom(ftr):
    for atom in Pharm.getSubstructure(ftr).atoms:
        if Chem.getType(atom) != Chem.AtomType.H:
            return atom

    return None
    
def postprocPharmacophore(ph4, is_lig):
    tot_hyd = 0.0
    
    for ftr in ph4:
        ftr_type = Pharm.getType(ftr)
        
        if ftr_type == Pharm.FeatureType.HYDROPHOBIC:
            hyd = Pharm.getHydrophobicity(ftr)
            tot_hyd += hyd
            
            Pharm.setWeight(ftr, hyd)
            
        elif ftr_type == Pharm.FeatureType.H_BOND_ACCEPTOR:
            acc_atom = getFeatureAtom(ftr)

            if is_lig:
                sybyl_type = Chem.getSybylType(acc_atom)

                if sybyl_type in SYBYL_TO_HBA_TYPE:
                    Pharm.setType(ftr, SYBYL_TO_HBA_TYPE[sybyl_type])
                else:
                    print('Unsupported sybyl H-acceptor type:', str(sybyl_type), file=sys.stderr)
                    
            else:
                atom_type = Chem.getType(acc_atom)

                if atom_type in ATOM_TO_HBA_TYPE:
                    Pharm.setType(ftr, ATOM_TO_HBA_TYPE[atom_type])
                    
        elif ftr_type == Pharm.FeatureType.H_BOND_DONOR:
            don_atom = getFeatureAtom(ftr)

            if is_lig:
                sybyl_type = Chem.getSybylType(don_atom)

                if sybyl_type in SYBYL_TO_HBD_TYPE:
                    Pharm.setType(ftr, SYBYL_TO_HBD_TYPE[sybyl_type])
                else:
                    print('Unsupported sybyl H-donor type:', str(sybyl_type), file=sys.stderr)

            else:
                atom_type = Chem.getType(don_atom)

                if atom_type in ATOM_TO_HBD_TYPE:
                    Pharm.setType(ftr, ATOM_TO_HBD_TYPE[atom_type])
                    
    return tot_hyd
    
def getHDonorFeatureCount(ph4, ftr_type):
    count = 0

    for ftr in ph4:
        if Pharm.getType(ftr) != ftr_type:
            continue

        don_atom = getFeatureAtom(ftr)
        count += MolProp.getAtomCount(don_atom, don_atom.getMolecule(), Chem.AtomType.H)
            
    return count

def getRotorBondCount(mol):
    count = 0

    for bond in mol.bonds:
        if ConfGen.isRotatableBond(bond, mol, True):
            count += 1

    return count
     
def processComplex(pdb_code, comp_data_dir, out_file, save_ph4s, exact, dyn_hbd):
    print('Processing complex %s...' % pdb_code)

    sdf_reader = Chem.FileSDFMoleculeReader(comp_data_dir + '/' + pdb_code + '_ligand.sdf')
    ligand = Chem.BasicMolecule()

    if not sdf_reader.read(ligand):
        print('!! While processing complex %s: reading ligand SD-file failed' % pdb_code, file=sys.stderr)
        return
        
    pdb_reader = Biomol.FilePDBMoleculeReader(comp_data_dir + '/' + pdb_code + '_protein.pdb')
    protein = Chem.BasicMolecule()

    if not pdb_reader.read(protein):
        print('!! While processing complex %s: reading protein PDB-file failed' % pdb_code, file=sys.stderr)
        return
    
    checkProtein(pdb_code, protein)
    
    lig_logp = prepareMolecule(ligand, True)
    prepareMolecule(protein, False)

    lig_env = Chem.Fragment()

    Biomol.extractEnvironmentResidues(ligand, protein, lig_env, Chem.Atom3DCoordinatesFunctor(), LIG_ENV_MAX_RADIUS, False)
    
    removeNonStdResidues(pdb_code, lig_env)

    Chem.perceiveSSSR(lig_env, True)
 
    lig_ph4 = Pharm.BasicPharmacophore()
    env_ph4 = Pharm.BasicPharmacophore()

    generatePharmacophore(ligand, lig_ph4, dyn_hbd or not exact, True)
    generatePharmacophore(lig_env, env_ph4, dyn_hbd, False)

    if save_ph4s:
        Pharm.FilePMLFeatureContainerWriter(comp_data_dir + '/' + pdb_code + '_lig_ph4.pml').write(lig_ph4)
        Pharm.FilePMLFeatureContainerWriter(comp_data_dir + '/' + pdb_code + '_env_ph4.pml').write(env_ph4)

    tot_lig_hyd = postprocPharmacophore(lig_ph4, True)
    postprocPharmacophore(env_ph4, False)
        
    line = pdb_code

    for ftr_type, ftr_name in FEATURE_NAMES.items():
        count = 0

        if (dyn_hbd or not exact) and ftr_type in LIG_HBD_FEATURE_TYPES:
            count = getHDonorFeatureCount(lig_ph4, ftr_type)
        else:
            if ftr_type in ENV_HBD_FEATURE_TYPES and ftr_type not in LIG_HBD_FEATURE_TYPES:
                continue

            if ftr_type in ENV_HBA_FEATURE_TYPES and ftr_type not in LIG_HBA_FEATURE_TYPES:
                continue
            
            count = Pharm.getFeatureCount(lig_ph4, ftr_type)
            
        line += ', ' + str(count)

    line += ', ' + str(tot_lig_hyd)
    line += ', ' + str(lig_logp)
#    line += ', ' + str(MolProp.calcLogS(ligand))
#    line += ', ' + str(getRotorBondCount(ligand))
            
    for ftr_type in ENV_HBA_FEATURE_TYPES:
        scores = calcEnvHDonorAcceptorOccupation(ligand, ftr_type, env_ph4, True)
        line += ', ' + str(scores[1])

    for ftr_type in ENV_HBD_FEATURE_TYPES:
        scores = calcEnvHDonorAcceptorOccupation(ligand, ftr_type, env_ph4, False)
        line += ', ' + str(scores[0])

    for ia_type in INTERACTION_TYPES:
        scores = calcScoresForInteraction(ia_type, lig_ph4, env_ph4, exact)

        if ia_type[3]:
            line += ', ' + str(scores[1])
        else:
            line += ', ' + str(scores[0])
             
    line += ', ' + str(calcElectrostaticInteractionEnergy(ligand, lig_env))
#    line += ', ' + str(calcLJVdWInteractionEnergy(ligand, lig_env))
    line += ', ' + str(calcMorseVdWInteractionEnergy(ligand, lig_env))

    out_file.write(line + '\n')
    out_file.flush()
    
def outputColNames(out_file):
    out_file.write('PDB code')

    for ftr_type, ftr_name in FEATURE_NAMES.items():
        if ftr_type in ENV_HBD_FEATURE_TYPES and ftr_type not in LIG_HBD_FEATURE_TYPES:
            continue

        if ftr_type in ENV_HBA_FEATURE_TYPES and ftr_type not in LIG_HBA_FEATURE_TYPES:
            continue

        out_file.write(', ' + ftr_name)

    out_file.write(', TOT_HYD')
    out_file.write(', LOGP')
#    out_file.write(', LOGS')
#    out_file.write(', ROT')

    for ftr_type in ENV_HBA_FEATURE_TYPES:
        #out_file.write(', ENV_' + FEATURE_NAMES[ftr_type] + '_OCC_SUM')
        out_file.write(', ENV_' + FEATURE_NAMES[ftr_type] + '_OCC_MAX')

    for ftr_type in ENV_HBD_FEATURE_TYPES:
        out_file.write(', ENV_' + FEATURE_NAMES[ftr_type] + '_OCC_SUM')
        #out_file.write(', ENV_' + FEATURE_NAMES[ftr_type] + '_OCC_MAX')

    for ia_type in INTERACTION_TYPES:
        ia_name = FEATURE_NAMES[ia_type[0]] + '-' + FEATURE_NAMES[ia_type[1]]

        if ia_type[3]:
            out_file.write(', ' + ia_name + '_MAX')
        else:
            out_file.write(', ' + ia_name + '_SUM')
    
    out_file.write(', ES, VDW\n')
    
def process(args):
    out_file = open(args.out_csv_file[0], 'w')

    outputColNames(out_file)
    
    for pdb_code in os.listdir(args.complex_data_dir[0]):
        comp_data_dir = os.path.join(args.complex_data_dir[0], pdb_code)

        if os.path.isfile(comp_data_dir): # sanity check
            continue

        try:
            df = processComplex(pdb_code, comp_data_dir, out_file, args.output_ph4s, args.exact, args.dyn_hbd)

        except Exception as e:
            print('!! Processing complex %s failed: ' % pdb_code, e, file=sys.stderr)

    out_file.close()
        
if __name__ == '__main__':
    process(parseArguments())
