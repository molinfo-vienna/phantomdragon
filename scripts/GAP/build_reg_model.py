# -*- mode: python; tab-width: 4 -*- 

## 
# build_reg_model.py 
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

import sys
import argparse

import CDPL.Math as Math


def parseArguments():
    parser = argparse.ArgumentParser(description='Builds a regression model for GRAIL-based affinity prediction.')
    
    parser.add_argument('-a',
                        dest='aff_values',
                        required=True,
                        help='[Required] The affinity value CSV-file.',
                        nargs=1)
    parser.add_argument('-i',
                        dest='aff_col_idx',
                        required=True,
                        help='[Required] The column index of the affinity value in the CSV-file.',
                        type=int,
                        nargs=1)
    parser.add_argument('-d',
                        dest='descr_file',
                        required=True,
                        help='[Required] The CSV-file containing the GRAIL-descriptors.',
                        nargs=1)
    parser.add_argument('-o',
                        dest='coeff_file',
                        required=False,
                        help='[Optional] The regression coefficients output file.',
                        nargs=1,
                        default=None)
    parser.add_argument('-x',
                        dest='excl_file',
                        required=False,
                        help='[Optional] A file containing PDB-codes (one per line) to exclude from model-building.',
                        nargs=1,
                        default=None)
    parser.add_argument('-p',
                        dest='poly',
                        required=False,
                        help='[Optional] Perform polynomial regression.',
                        action='store_true',
                        default=False)
 
    return parser.parse_args()

def process(args):
    print('Calculation regression coefficients...')

    aff_value_csv = open(args.aff_values[0], 'r')
    aff_value_idx = args.aff_col_idx[0]
    descrs_csv = open(args.descr_file[0], 'r')
    descr_map = {}
    
    ftr_names = descrs_csv.readline().split(',')
    num_ftrs = len(ftr_names)

    for i in range(1, num_ftrs):
        ftr_names[i] = ftr_names[i].strip()

        if args.poly:
            ftr_names.append(ftr_names[i] + '_SQRD')

    ftr_names[0] = 'Offset'

    while True:
        line = descrs_csv.readline()

        if not line:
            break

        cols = line.split(',')
        pdb_code = cols[0].strip()

        descr_map[pdb_code] = cols

    mlrm = Math.DMLRModel()
    descr_vec = Math.DVector()

    excl_list = set()

    if args.excl_file:
        for pdb_code in open(args.excl_file[0], 'r').readlines():
            excl_list.add(pdb_code.strip())
    
    while True:
        line = aff_value_csv.readline()

        if not line:
            break

        cols = line.split(' ')
        pdb_code = cols[0].strip()

        if pdb_code in excl_list:
            continue
        
        if pdb_code in descr_map:
            aff_value = float(cols[1 + aff_value_idx])
            descr_line = descr_map[pdb_code]

            if args.poly:
                descr_vec.resize(len(descr_line) * 2 - 1, 0.0)
            else:
                descr_vec.resize(len(descr_line))

            descr_vec[0] = 1.0

            for i in range(1, len(descr_line)):
                descr_vec[i] = float(descr_line[i].strip())

                if args.poly:
                    descr_vec[i + len(descr_line) - 1] = descr_vec[i] ** 2

            mlrm.addXYData(descr_vec, aff_value)

    mlrm.buildModel()
    mlrm.calcStatistics()
    
    print('Num. data points =', mlrm.yValues.size)
    print('Chi squared =', mlrm.chiSquare)
    print('Goodness of fit =', mlrm.goodnessOfFit)
    print('Correlation coefficient =', mlrm.correlationCoefficient)
    print('Standard deviation =', mlrm.standardDeviation)
    print('Num. descr. features =', descr_vec.getSize())

    if args.coeff_file:
        coeffs_file = open(args.coeff_file[0], 'w')
        new_line = True

        for i in range(0, len(mlrm.coefficients)):
            coeff = str(mlrm.coefficients[i])

            if not new_line:
                coeffs_file.write(', ')

            coeffs_file.write(coeff)
            new_line = False

            if (i + 1) % 5 == 0:
                coeffs_file.write(',\n')
                new_line = True

            #print(ftr_names[i] + ':\t' + coeff)

    print('Done!')

if __name__ == '__main__':
    process(parseArguments())
