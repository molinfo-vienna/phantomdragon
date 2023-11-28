# -*- mode: python; tab-width: 4 -*- 

## 
# pred_aff_pdb_ligands.py 
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

from os import path

import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math
import CDPL.GRAIL as GRAIL


LIG_ENV_MAX_RADIUS = 21.0
REMOVE_NON_STD_RESIDUES = True


def parseArguments():
    parser = argparse.ArgumentParser(description='Performs GRAIL affinity predictions for a PDB-file and set of input ligands.')
    parser.add_argument('-p',
                        dest='pdb_file',
                        required=True,
                        help='[Required] The receptor PDB-file.',
                        nargs=1)
    parser.add_argument('-l',
                        dest='lig_file',
                        required=True,
                        help='[Required] The file providing the ligands.',
                        nargs=1)
    parser.add_argument('-o',
                        dest='out_csv_file',
                        required=True,
                        help='[Required] The path of the output CSV-file containing the predicted affinity values calculated for each input ligand.',
                        nargs=1)
 
    return parser.parse_args()

def removeNonStdResidues(pdb_file, protein):
    residues = Biomol.ResidueList(protein)
    
    for res in residues:
        is_std_res = Biomol.ResidueDictionary.isStdResidue(Biomol.getResidueCode(res))
        
        if is_std_res and res.numAtoms < 5:
            print('!! While processing %s: isolated standard residue fragment of size %s found' % (pdb_file, str(res.numAtoms)), file=sys.stderr)
            protein -= res
            
        elif REMOVE_NON_STD_RESIDUES and not is_std_res:
            if res.numAtoms == 1 and Chem.AtomDictionary.isMetal(Chem.getType(res.atoms[0])):
                continue

            protein -= res

    Chem.clearSSSR(protein)

def checkProtein(pdb_file, protein):
    for atom in protein.atoms:
        if Chem.getType(atom) == Chem.AtomType.H and atom.numAtoms == 0:
            print('!! While processing %s: isolated hydrogen atom encountered' % pdb_file, file=sys.stderr)
             
        elif Chem.getType(atom) == Chem.AtomType.UNKNOWN:
            print('!! While processing %s: atom of unknown element encountered' % pdb_file, file=sys.stderr)

def loadPDBFile(pdb_file):
    print('Loading PDB-file %s...' % path.basename(pdb_file))
    
    pdb_reader = Biomol.FilePDBMoleculeReader(pdb_file)
    protein = Chem.BasicMolecule()

    if not pdb_reader.read(protein):
        sys.exit('!! Reading PDB-file %s failed' % path.basename(pdb_file), file=sys.stderr)
    
    checkProtein(path.basename(pdb_file), protein)
    removeNonStdResidues(path.basename(pdb_file), protein)
    
    GRAIL.prepareForGRAILDescriptorCalculation(protein, True)

    return protein

def processComplex(protein, ligand, out_file, descr_calc, lig_idx):
    GRAIL.prepareForGRAILDescriptorCalculation(ligand, True)

    lig_env = Chem.Fragment()

    Biomol.extractEnvironmentResidues(ligand, protein, lig_env, Chem.Atom3DCoordinatesFunctor(), LIG_ENV_MAX_RADIUS, False)
    Chem.extractSSSRSubset(protein, lig_env, True)

    line = Chem.getName(ligand)

    if not line:
        line = str(lig_idx)

    descr = Math.DVector()
    lig_atom_coords = Math.Vector3DArray()

    Chem.get3DCoordinates(ligand, lig_atom_coords)

    descr_calc.initTargetData(lig_env, Chem.Atom3DCoordinatesFunctor())
    descr_calc.initLigandData(ligand)

    descr_calc.calculate(lig_atom_coords, descr)
 
    aff_calc = GRAIL.BindingAffinityCalculator()

    line += ', ' + str(aff_calc(descr, GRAIL.BindingAffinityCalculator.PKI))
    line += ', ' + str(aff_calc(descr, GRAIL.BindingAffinityCalculator.PKD))
    line += ', ' + str(aff_calc(descr, GRAIL.BindingAffinityCalculator.PKD_PKI))

    out_file.write(line + '\n')
    out_file.flush()

def outputColNames(out_file):
    out_file.write('Ligand, pKi, pKd, pKi/pKd\n')
    out_file.flush()
    
def process(args):
    protein = loadPDBFile(args.pdb_file[0])
    
    lig_reader = Chem.MoleculeReader(args.lig_file[0])

    Chem.setMultiConfImportParameter(lig_reader, False)

    ligand = Chem.BasicMolecule()

    out_file = open(args.out_csv_file[0], 'w')
    descr_calc = GRAIL.GRAILDescriptorCalculator()
    i = 1
    
    outputColNames(out_file)    

    print('Calculating affinities for ligands in %s...' % path.basename(args.lig_file[0]))
     
    while lig_reader.read(ligand):
        try:
            processComplex(protein, ligand, out_file, descr_calc, i)
            i += 1

        except Exception as e:
            print('!! Processing complex failed: ', e, file=sys.stderr)

    out_file.close()

    print('Done!')
    
if __name__ == '__main__':
    process(parseArguments())
