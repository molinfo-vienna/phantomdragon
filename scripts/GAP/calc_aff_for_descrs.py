# -*- mode: python; tab-width: 4 -*- 

## 
# calc_aff_for_descrs.py 
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
import numpy


def parseArguments():
    parser = argparse.ArgumentParser(description='Calculates GRAIL affinity values for a list of GRAIL descriptors.')
    parser.add_argument('-d',
                        dest='descr_file',
                        required=True,
                        help='[Required] The CSV-file containing the GRAIL descriptors.',
                        nargs=1)
    parser.add_argument('-c',
                        dest='coeff_file',
                        required=True,
                        help='[Required] The file with the affinity model regression coefficients.',
                        nargs=1) 
    parser.add_argument('-t',
                        dest='aff_type',
                        required=True,
                        help='[Required] The type of affinity value to calculate.',
                        nargs=1) 
    parser.add_argument('-o',
                        dest='out_csv_file',
                        required=True,
                        help='[Required] The path of the output CSV-file containing the affinities calculated for each GRAIL descriptor.',
                        nargs=1)

    return parser.parse_args()

def readCoefficients(coeff_file):
    coeff_list = []

    for line in coeff_file.readlines():
        values = line.split(',')

        for value in values:
            value = value.strip()

            if value:
                coeff_list.append(float(value.strip()))

    return numpy.array(coeff_list)
    
def process(args):
    print('Calculating affinity values...')

    descr_file = open(args.descr_file[0])
    coeff_file = open(args.coeff_file[0])
    aff_type = args.aff_type[0]
    out_file = open(args.out_csv_file[0], 'w')
   
    coeffs = readCoefficients(coeff_file)

    out_file.write(descr_file.readline().split(',')[0] + ',' + aff_type + '\n')
    
    while True:
        descr_line = descr_file.readline()

        if not descr_line:
            break

        descr_values = descr_line.split(',')
        
        out_file.write(descr_values[0] + ',')

        descr_values[0] = 1.0

        for i in range(1, len(descr_values)):
            descr_values[i] = float(descr_values[i].strip())

        if len(coeffs) > len(descr_values):
            for i in range(1, len(descr_values)):
                descr_values.append(descr_values[i] ** 2)

        aff = numpy.inner(coeffs, descr_values)

        out_file.write(str(aff) + '\n')

    cols = descr_line

    out_file.close()

    print('Done!')
    
if __name__ == '__main__':
    process(parseArguments())
