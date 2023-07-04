#!/usr/bin/env python

# Copyright 2022 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import os

import numpy as np

from rdkit import Chem
from rdkit.Chem.Descriptors import HeavyAtomMolWt
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import connected_components
import itertools

from dm_job_utilities.dm_log import DmLog
import utils


def _powerset(iterElems, min_num_elements=2, max_num_emements=None, include_full=False, combinations_instead_permutations=True):
    s = list(iterElems)
    if max_num_emements is None:
        last = len(s)+1 if include_full else len(s)
    else:
        last = max_num_emements+1
    func = itertools.combinations if combinations_instead_permutations else itertools.permutations
    return itertools.chain.from_iterable(func(s, r) for r in range(min_num_elements, last))


def get_permutations(list_of_fragments: Chem.Mol, min_num=2, max_num=3, min_dist=1.5):

    combinations_iter = _powerset(list_of_fragments, min_num_elements=min_num, max_num_emements=max_num,
                                  include_full=False, combinations_instead_permutations=True)
    for frag_group in combinations_iter:
        coords = [mol.GetConformer().GetPositions() for mol in frag_group]
        adj_mat = np.zeros((len(coords), len(coords)))
        for i, coord0 in enumerate(coords):
            for j, coord1 in enumerate(coords):
                dists = cdist(coord0, coord1)
                if dists.min() < min_dist:
                    adj_mat[i, j] = adj_mat[j, i] = 1
        n_components, labels = connected_components(csgraph=adj_mat, directed=False)
        if n_components == 1:
            yield frag_group


def generate_combinations(inputs, output_dir, min_dist, min_num, max_num):

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    else:
        output_dir = './'

    mols = utils.read_molecules(inputs)
    DmLog.emit_event('Found', len(mols), 'molecules')

    count = 0
    for perm in get_permutations(mols, min_dist=min_dist, min_num=min_num, max_num=max_num):
        count += 1
        indexes = [str(mols.index(m)) for m in perm]
        DmLog.emit_event('Writing', indexes)
        fname = os.path.join(output_dir, '_'.join(indexes) + '.sdf')
        with Chem.SDWriter(fname) as writer:
            for m in perm:
                writer.write(m)

    DmLog.emit_event('Generated', count, 'combinations')


def main():

    # Example:
    #   ./prep_compatible_frags.py -i data/x0032_0A.mol data/x0103_0A.mol -o output
    # or
    #   ./prep_compatible_frags.py -i data/mpro_hits.sdf -o output

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare compatible fragments')
    parser.add_argument('-i', '--inputs', nargs='+', required=True, help="Molfiles or SD-files with fragments")
    parser.add_argument('-o', '--output-dir', help="Location of the output dir. If not defined then ./ is used")
    parser.add_argument('--min-dist', type=float, default=1.5, help="Minimum distance threshold")
    parser.add_argument('--max-num', type=int, default=2, help="Maximum permutation size")
    parser.add_argument('--min-num', type=int, default=2, help="Minimum permutation size")

    args = parser.parse_args()
    DmLog.emit_event("prep_compatible_frags: ", args)

    generate_combinations(args.inputs, args.output_dir, args.min_dist, args.min_num, args.max_num)


if __name__ == "__main__":
    main()
