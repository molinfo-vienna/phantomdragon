## 
# cdf_generation_md.py 
#
# Copyright (C) 2023 Thomas A. Seidel, Christian Fellinger and Arthur Garon
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
from collections import defaultdict

import MDAnalysis

import CDPL.Base as Base
import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Math as Math


def loadCDFMolecule(fname):
    mol = Chem.BasicMolecule()
    cdf_reader = Chem.FileCDFMoleculeReader(fname)

    try:
        if not cdf_reader.read(mol):
            sys.exit('!! Could not load CDF file: ' + fname)
    except:
        sys.exit('!! Could not load CDF file: ' + fname)
            
    Chem.calcImplicitHydrogenCounts(mol, False)
    Chem.perceiveHybridizationStates(mol, False)
    Chem.setAtomSymbolsFromTypes(mol, False)
    Chem.perceiveSSSR(mol, False)
    Chem.setRingFlags(mol, False)
    Chem.setAromaticityFlags(mol, False)

    return mol

def mergeCDFMolecule(fname, chunk_number):
    print('Merging %s CDF-file chunks...' % str(chunk_number))
    
    main_cdf = loadCDFMolecule(fname)

    for chunk_index in range(1, int(chunk_number)):
        tmp_cdf = loadCDFMolecule(os.path.dirname(fname) +'/'+ os.path.basename(fname).split('_chunk_')[0] + '_chunk_' + str(chunk_index) + '.cdf')
        coords_func = Chem.AtomConformer3DCoordinatesFunctor(0)

        for atom_index in range(tmp_cdf.numAtoms):
            main_coords_array = Chem.get3DCoordinatesArray(main_cdf.getAtom(atom_index))
            tmp_coords_array = Chem.get3DCoordinatesArray(tmp_cdf.getAtom(atom_index))

            for position_index in range(tmp_coords_array.size):
                main_coords_array.addElement(tmp_coords_array[position_index])

    main_file = fname.split('_chunk_')[0] + '.cdf'
    cdf_writer = Chem.FileCDFMolecularGraphWriter(main_file)

    try:
        if not cdf_writer.write(main_cdf):
            sys.exit('!! Could not write merged CDF file: ' + main_file)
    except:
        sys.exit('!! Could not write merged CDF file: ' + main_file)
        
    for chunk_index in range(chunk_number):
        os.remove(os.path.dirname(fname) +'/'+ os.path.basename(fname).split('_chunk_')[0] + '_chunk_' + str(chunk_index) + '.cdf')

    cdf_writer.close()

    print('Done!')

def generateCDFMolecule(top_file, traj_files, output, name, chunk_size=500):
    u = MDAnalysis.Universe(top_file, traj_files)
    chunks_start = [0]
    chunks_end = []

    print('Loaded trajectory with %s frames...' % str(len(u.trajectory)));
    
    if chunk_size == 0 or chunk_size > len(u.trajectory):
        chunk_size = len(u.trajectory)

    for i in range(len(u.trajectory)):
        if i % chunk_size == 0 and i != 0 and i != (len(u.trajectory) - 1):
            chunks_end.append(i)
            chunks_start.append(i)
        elif i == (len(u.trajectory) - 1):
            chunks_end.append(i + 1)

    cdf_mol = Chem.BasicMolecule()
    waters = {}
    
    if top_file.endswith('.pdb'):
        pdb_str = open(top_file, 'r').read().replace('WAT', 'HOH').replace('HIE', 'HIS').replace('SPC X', 'HOH X')
        pdb_reader = Biomol.PDBMoleculeReader(Base.StringIOStream(pdb_str))

        Biomol.setPDBApplyDictAtomBondingToNonStdResiduesParameter(pdb_reader, True)

        try:
            if not pdb_reader.read(cdf_mol):
                sys.exit('!! Could not read PDB file: ' + top_file)
        except:
            sys.exit('!! Could not read PDB file: ' + top_file)

        for atom in cdf_mol.atoms:
            Chem.setImplicitHydrogenCount(atom, 0)
            
            if Biomol.getResidueCode(atom) == 'HOH':
                res_id = Biomol.getResidueSequenceNumber(atom)
                 
                if res_id in waters:
                    waters[res_id].append(int(res_id))
                else:
                    waters[res_id] = [int(res_id)]
                    
            array = Math.Vector3DArray()
            array.resize(chunk_size, Math.Vector3D())

            Chem.set3DCoordinatesArray(atom, array)
    else:
        cdf_mol.reserveMemoryForAtoms(len(u.atoms))
        cdf_mol.reserveMemoryForBonds(len(u.bonds))

        # construct atoms
        for md_atom in u.atoms:
            atom = cdf_mol.addAtom()

            if md_atom.type == 'CL':
                Chem.setSymbol(atom, 'Cl')
            else:
                Chem.setSymbol(atom, md_atom.name[0])

            Biomol.setChainID(atom, md_atom.segid)

            if md_atom.resname == 'WAT' or md_atom.resname == 'TIP3':
                Biomol.setResidueCode(atom, 'HOH')
            elif md_atom.resname == 'HSD' or md_atom.resname == 'HIE':
                Biomol.setResidueCode(atom, 'HIS')
            else:
                Biomol.setResidueCode(atom, md_atom.resname)

            Biomol.setResidueSequenceNumber(atom, int(md_atom.resid))
            Biomol.setResidueAtomName(atom, md_atom.name)
            Biomol.setSerialNumber(atom, int(md_atom.id))

            if Biomol.getResidueCode(atom) == 'HOH':
                if md_atom.resid in waters:
                    waters[md_atom.resid].append(int(md_atom.id))
                else:
                    waters[md_atom.resid] = [int(md_atom.id)]

            # fix positive charge on arginin nitrogen
            if md_atom.resname == 'ARG' and md_atom.name == 'NH2':
                Chem.setFormalCharge(atom, 1)

            Chem.set3DCoordinates(cdf_mol.getAtom(int(md_atom.id)), md_atom.position.tolist())
            coords_array = Math.Vector3DArray()
            coords_array.resize(chunk_size, Math.Vector3D())

            Chem.set3DCoordinatesArray(atom, coords_array)

        Chem.setAtomTypesFromSymbols(cdf_mol, True)

        # construct bonds
        for md_bond in u.bonds:
            if not Chem.getType(cdf_mol.atoms[int(md_bond.atoms[0].index)]) == Chem.AtomType.H == Chem.getType(
                    cdf_mol.atoms[int(md_bond.atoms[1].index)]):
                cdf_mol.addBond(int(md_bond.atoms[0].index), int(md_bond.atoms[1].index))
   
        # make sane biomolecule
        for a in cdf_mol.atoms:
            Chem.setImplicitHydrogenCount(a, 0)

    for water in waters.values():
        if len(water) < 2:
            continue

        for atom_idx in water:
            if Chem.getType(cdf_mol.atoms[atom_idx]) == Chem.AtomType.O:
                if cdf_mol.atoms[atom_idx].numBonds > 1:
                    break

                for atom_idx2 in water:
                    if Chem.getType(cdf_mol.atoms[atom_idx2]) == Chem.AtomType.H:
                        cdf_mol.addBond(atom_idx, atom_idx2)

                break

    Chem.perceiveSSSR(cdf_mol, True)
    Chem.setRingFlags(cdf_mol, True)
    Chem.perceiveBondOrders(cdf_mol, True)
    Chem.perceiveHybridizationStates(cdf_mol, True)
    Chem.setAromaticityFlags(cdf_mol, True)
    Chem.calcFormalCharges(cdf_mol, True)

    # read timesteps
    u = MDAnalysis.Universe(top_file, traj_files)

    for i in range(len(chunks_start)):
        print('Processing frames %s-%s...' % (str(chunks_start[i]), str(chunks_end[i] - 1)))
        
        tmp_array = defaultdict(lambda: [])
        for ts in u.trajectory[chunks_start[i]: chunks_end[i]]:
            for md_atom in u.atoms:
                tmp_array[int(md_atom.index)].append([x for x in md_atom.position])

        for index in list(tmp_array.keys()):
            coords_array = Chem.get3DCoordinatesArray(cdf_mol.getAtom(index))
            for element_id, element in enumerate(tmp_array[index]):
                for dim in range(3):
                    coords_array[element_id][dim] = float(element[dim])

        if (i+1)*chunk_size > len(u.trajectory):
            over_size = (i+1)*chunk_size - len(u.trajectory)

            for md_atom in cdf_mol.atoms:
                coords_array = Chem.get3DCoordinatesArray(md_atom)
                for y in range(over_size):
                    coords_array.popLastElement()

        tmp_output = output + name + '_chunk_' + str(i) + '.cdf'

        try:
            if not Chem.FileCDFMolecularGraphWriter(tmp_output).write(cdf_mol):
                sys.exit('!! Could not write CDF file: ' + tmp_output)
        except:
            sys.exit('!! Could not write CDF file: ' + tmp_output)

    print('Done!')
            
    return len(chunks_start)

def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate a cdf file (CDPKit format) from a topology and one or more trajectory file(s).')
    parser.add_argument('-trj',
                        dest='trajectory',
                        required=True,
                        help='[Required] The path of the trajectory file, can also handle multiple trajectories of the same system.',
                        nargs="*")
    parser.add_argument('-top',
                        dest='topology',
                        required=True,
                        help='[Required] The path of the topology file.',
                        nargs=1)
    parser.add_argument('-n',
                        dest='cdf_name',
                        help='[Optional] The name of the cdf file (Default: name of the trajectory file).',
                        nargs=1,
                        default=None)
    parser.add_argument('-cs',
                        dest='chunk_size',
                        help='[Optional] The number of frames to consider per chunk (default: 500).',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the cdf file will be generated (default: current directory).',
                        nargs=1,
                        default=None)

    return parser.parse_args()

if __name__ == '__main__':
    args = parseArguments()

    trajectory = args.trajectory
    topology = args.topology[0]

    if args.cdf_name is None:
        cdf_name = os.path.basename(trajectory[0])[:-4]
    else:
        cdf_name = args.cdf_name[0]

    if args.chunk_size is None:
        chunk_size = 500
    else:
        chunk_size = int(args.chunk_size[0])

    if args.output is None:
        output = './'
    else:
        output = args.output[0] + '/'

    if output[-1] != '/':
        output += '/'

    chunk_number = generateCDFMolecule(topology, trajectory, output, cdf_name, chunk_size=chunk_size)
    mergeCDFMolecule(output + cdf_name + '_chunk_0.cdf', chunk_number)
