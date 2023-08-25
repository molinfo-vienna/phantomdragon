# -*- mode: python; tab-width: 4 -*- 

## 
# calc_pred_stats.py 
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
import math


def parseArguments():
    parser = argparse.ArgumentParser(description='Calculates GRAIL affinity value prediction statistics.')
    parser.add_argument('-p',
                        dest='pred_file',
                        required=True,
                        help='[Required] The CSV-file containing the predicted affinities.',
                        nargs=1)
    parser.add_argument('-r',
                        dest='ref_file',
                        required=True,
                        help='[Required] The CSV-file containing the reference affinities.',
                        nargs=1)
    parser.add_argument('-i',
                        dest='ref_col_idx',
                        required=True,
                        help='[Required] The column index of the affinity value in the reference CSV-file.',
                        type=int,
                        nargs=1)

    return parser.parse_args()

def readAffinities(in_file, is_ref_file = False, aff_col_idx = 1):
    values = {}
    first_line = True
    
    while True:
        line = in_file.readline()

        if not line:
            break

        if first_line and not is_ref_file:
            first_line = False
            continue

        if is_ref_file:
            tokens = line.split(' ')
            values[tokens[0].strip()] = float(tokens[aff_col_idx].strip())
        else:
            tokens = line.split(',')
            values[tokens[0].strip()] = float(tokens[1].strip())
        
    return values
        
def process(args):
    print('Calculating statistics...')

    pred_affs = readAffinities(open(args.pred_file[0]))
    ref_affs = readAffinities(open(args.ref_file[0]), True, args.ref_col_idx[0])
    pred_ref_aff_pairs = []

    for pdb_code, ref_aff in ref_affs.items():
        if pdb_code in pred_affs:
            pred_ref_aff_pairs.append((pred_affs[pdb_code], ref_aff))

    pred_mean = 0.0
    ref_mean = 0.0
    std_err = 0.0

    for p, r in pred_ref_aff_pairs:
        pred_mean += p
        ref_mean += r
        std_err += (p - r) ** 2

    pred_mean /= len(pred_ref_aff_pairs)
    ref_mean /= len(pred_ref_aff_pairs)

    sxx = 0
    syy = 0
    sxy = 0

    for p, r in pred_ref_aff_pairs:
        xt = r - ref_mean
        yt = p - pred_mean

        sxx += xt * xt
        syy += yt * yt
        sxy += xt * yt

    r = sxy / math.sqrt(sxx * syy)
    std_err = math.sqrt(std_err / (len(pred_ref_aff_pairs)))

    print('Num. points: ' + str(len(pred_ref_aff_pairs)))
    print('Standard error: ' + str(std_err))
    print('Pearson r: ' + str(r))
    
    print('Done!')
    
if __name__ == '__main__':
    process(parseArguments())
