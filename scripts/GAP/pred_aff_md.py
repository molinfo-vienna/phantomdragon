# -*- mode: python; tab-width: 4 -*- 

## 
# pred_aff_md.py 
#
# Copyright (C) 2023 Thomas A. Seidel <thomas.seidel@univie.ac.at>
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

import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math
import CDPL.GRAIL as GRAIL


LIG_ENV_MAX_RADIUS = 25.0
REMOVE_NON_STD_RESIDUES = True
NON_STD_RESIDUE_MAX_ATOM_COUNT = 8


def parseArguments():
    parser = argparse.ArgumentParser(description='Calculates GRAIL affinity prediction descriptors for MD-generated ligand/receptor conformations.')
    parser.add_argument('-i',
                        dest='traj_file',
                        required=True,
                        help='[Required] The CDF-file generated from the MD trajectory.',
                        nargs=1)
    parser.add_argument('-l',
                        dest='lig_tlc',
                        required=True,
                        help='[Required] The three-letter code of the ligand.',
                        nargs=1) 
    parser.add_argument('-o',
                        dest='out_csv_file',
                        required=True,
                        help='[Required] The path of the output CSV-file containing the descriptor values calculated for each MD frame.',
                        nargs=1)
    parser.add_argument('-c',
                        dest='norm_chgs',
                        help='[Optional] Change protonation of acidic/basic groups to a state likely at pH7 (default: false)·',
                        action='store_true',
                        default=False)
    return parser.parse_args()

def removeNonStdResidues(lig_env):
    residues = Biomol.ResidueList(lig_env)
    
    for res in residues:
        is_std_res = Biomol.ResidueDictionary.isStdResidue(Biomol.getResidueCode(res))
        
        if is_std_res and res.numAtoms < 5:
            print('!! Isolated standard residue fragment of size %s found' % str(res.numAtoms), file=sys.stderr)
            lig_env -= res
            
        elif REMOVE_NON_STD_RESIDUES and not is_std_res and res.numAtoms <= NON_STD_RESIDUE_MAX_ATOM_COUNT:
            if res.numAtoms == 1 and Chem.AtomDictionary.isMetal(Chem.getType(res.atoms[0])):
                continue

            lig_env -= res

def checkComplex(comp):
    for atom in comp.atoms:
        if Chem.getType(atom) == Chem.AtomType.H and atom.numAtoms == 0:
            print('!! Isolated hydrogen atom encountered', file=sys.stderr)
             
        elif Chem.getType(atom) == Chem.AtomType.UNKNOWN:
            print('!! Atom of unknown element encountered', file=sys.stderr)

def processTrajectory(traj_file, lig_tlc, norm_chgs, out_file):
    print('Processing trajectory file %s...' % traj_file)

    comp = Chem.BasicMolecule()
    cdf_reader = Chem.FileCDFMoleculeReader(traj_file)

    try:
        if not cdf_reader.read(comp):
            sys.exit('!! Could not load trajectory file ' + traj_file)
    except Exception as e:
        sys.exit('!! Could not load trajectory file {}: {}' % (traj_file, str(e)))
    
    checkComplex(comp)
    removeNonStdResidues(comp)
    GRAIL.prepareForGRAILDescriptorCalculation(comp, norm_chgs)
    
    ligand = Chem.Fragment()
    lig_env = Chem.Fragment()
    
    for atom in comp.atoms:
        if Biomol.getResidueCode(atom) == lig_tlc:
            Biomol.extractResidueSubstructure(atom, comp, ligand, False)
            break

    if ligand.numAtoms == 0:
        sys.exit('!! Could not find ligand {}'.format(lig_tlc))
        
    Biomol.extractEnvironmentResidues(ligand, comp, lig_env, Chem.Atom3DCoordinatesFunctor(), LIG_ENV_MAX_RADIUS, False)
    Chem.extractSSSRSubset(comp, lig_env, True)
    Chem.extractSSSRSubset(comp, ligand, True)

    descr_calc = GRAIL.GRAILDescriptorCalculator()
    aff_calc = GRAIL.BindingAffinityCalculator()
    
    descr_calc.initLigandData(ligand)
    
    descr = Math.DVector()
    lig_atom_coords = Math.Vector3DArray()
    num_confs = Chem.getNumConformations(comp)

    out_file.write('Frame, pKi, pKd, pKi/pKd\n')
    out_file.flush();
    
    print('Calculating ligand binding affinities for %s frames...' % str(num_confs))

    first_call = True
    
    for i in range(0, num_confs):
        line = str(i)
    
        Chem.getConformation(ligand, i, lig_atom_coords)

        descr_calc.initTargetData(lig_env, Chem.AtomConformer3DCoordinatesFunctor(i), first_call)
        descr_calc.calculate(lig_atom_coords, descr, first_call)

        line += ', ' + str(aff_calc(descr, GRAIL.BindingAffinityCalculator.PKI))
        line += ', ' + str(aff_calc(descr, GRAIL.BindingAffinityCalculator.PKD))
        line += ', ' + str(aff_calc(descr, GRAIL.BindingAffinityCalculator.PKD_PKI))
        
        out_file.write(line + '\n')
        out_file.flush()

        first_call = False
    
def process(args):
    out_file = open(args.out_csv_file[0], 'w')
    
    try:
        df = processTrajectory(args.traj_file[0], args.lig_tlc[0], args.norm_chgs, out_file)

    except Exception as e:
        print('!! Processing trajectory file %s failed: ' % args.traj_file[0], e, file=sys.stderr)

    out_file.close()

    print('Done!')
    
if __name__ == '__main__':
    process(parseArguments())
