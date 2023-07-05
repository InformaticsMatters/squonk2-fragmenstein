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

# TODO:
# - better handling of fragmenstein stdout/stderr
# - run Victor as subprocess?

from dm_job_utilities.dm_log import DmLog
from fragmenstein import Wictor

from rdkit import Chem
import gzip, argparse, time, traceback, os, random, string
import utils


job_id = ''.join(random.choice(string.ascii_lowercase) for i in range(16))


def create_victor(hits, pdb, work_dir, seed):
    """
    Create the Victor instance. It might be possible to re-use this for multiple placements, but the time taken to
    create the instance is negligible compared to the placement process.

    :param hits: The hits to align to (RDKit molecule)
    :param pdb: Path of the PDB file for the protein
    :param seed: Random seed for Monster. If not defined then a random seed will be generated
    :return: The Victor instance
    """
    if not seed:
        seed = random.randint(0, int(1e6))
    victor = Wictor(hits=hits, pdb_filename=pdb, monster_random_seed=seed)
    if work_dir:
        utils.log('Setting work_path to', work_dir)
        victor.work_path = work_dir
    return victor


def prepare_mol_for_writing(mol, index, id, protein, frag_mols, ddG, rmsd, removeHydrogens=True,
                            idxProp='IDX', ddgProp='DDG', rmsdProp='RMSD', smilesProp=None,
                            refMolsProp=None, fragIdField=None, proteinProp=None, proteinPropValue=None):

    mol_noh = Chem.RemoveHs(mol)
    if removeHydrogens:
        mol = mol_noh

    mol.SetProp('_Name', job_id + '-' + str(index))
    mol.SetProp(idxProp, str(id))
    mol.SetDoubleProp(ddgProp, ddG)
    mol.SetDoubleProp(rmsdProp, rmsd)
    if smilesProp:
        mol.SetProp(smilesProp, Chem.MolToSmiles(mol_noh))
    if refMolsProp and fragIdField:
        add_input_ids(mol, frag_mols, fragIdField, refMolsProp)
    if proteinProp:
        if proteinPropValue:
            mol.SetProp(proteinProp, proteinPropValue)
        else:
            mol.SetProp(proteinProp, protein)
    return mol


def add_input_ids(mol, inputs, nameProp, refMolsProp):
    refs = []
    for i, input in enumerate(inputs):
        if input.HasProp(nameProp):
            val = input.GetProp(nameProp)
            refs.append(val)
    mol.SetProp(refMolsProp, ','.join(refs))


def place_smiles(fragment_mols, smiles, protein, work_dir=None, name='fragmenstein', seed=None):
    """
    Place an individual smiles

    :param fragment_mols: The hits to align to (RDKit molecules)
    :param smiles: The SMILES to place
    :param protein: Path of the PDB file for the protein
    :param name: Job name
    :param seed: Random seed for Monster. If not defined then a random seed will be generated
    :return: The embedded molecules (RDKit molecule)
    """
    DmLog.emit_event('Placing SMILES', smiles)

    v = create_victor(fragment_mols, protein, work_dir, seed)
    v.place(smiles, long_name=name)
    mol = v.minimized_mol
    return mol, v.ddG, v.mrmsd.mrmsd


def read_smiles_file(smiles_file, readHeader=False, idColumn=1, delimiter=None):
    smiles = []
    ids = []
    if smiles_file.endswith('.gz'):
        reader = gzip.open(smiles_file, 'rt')
    else:
        reader = open(smiles_file, 'rt')

    if readHeader:
        header = reader.readline()

    mol_idx = 0
    while True:
        line = reader.readline()
        if not line:
            break
        else:
            mol_idx += 1
            line = line.strip()
            if delimiter is None:
                tokens = line.split()
            else:
                tokens = line.split(delimiter)
            stripped = []
            for token in tokens:
                stripped.append(token.strip())
            smiles.append(tokens[0])
            if len(tokens) >= idColumn:
                ids.append(tokens[idColumn])
            else:
                ids.append(str(mol_idx))

    return smiles, ids


def place_multiple_smiles(fragments, smiles, ids, protein, outputFileName, work_dir=None,
                    fragIdField=None, refMolsProp=None, removeHydrogens=False,
                    idxProp='IDX', ddgProp='DDG', rmsdProp='RMSD', smilesProp=None,
                    proteinProp=None, proteinPropValue=None, num_runs=1):
    """
    Place a list of SMILES, writing the embedded molecules to a SD file

    :param fragments: The fragments to align to (RDKit molecules)
    :param smiles: List of SMILES strings
    :param ids: Optional IDs for the SMILES. If not defined then the index is used as the ID
    :param protein: Path of the PDB file for the protein
    :param outputFileName: Name of the output SD file
    :param work_dir: directory where fragmenstein does its work
    :param fragIdField: field name for the ID of the input
    :param refMolsProp: Field name for the IDs of the input molecules
    :param idxProp: Name of the property in the output file with the smiles file index (default IDX)
    :param num_runs: Number of instantiations for each SMILES (default 1)
    :return: num_mols, num_embeddings, num_failures
    """

    if not os.path.isfile(protein):
        raise ValueError("Protein file not found: " + protein)

    frag_mols = utils.read_molecules(fragments)
    DmLog.emit_event('Read', len(frag_mols), 'fragments')

    utils.expand_path(outputFileName)

    count = 0
    num_placements = 0
    num_errors = 0
    mol_idx = 0
    with Chem.SDWriter(outputFileName) as writer:
        for smi in smiles:
            count += 1
            try:
                if ids and len(ids) >= mol_idx:
                    id = ids[mol_idx]
                else:
                    id = str(mol_idx + 1)
                DmLog.emit_event('Processing SMILES', id)
                for i in range(num_runs):
                    DmLog.emit_event('Running placement', i + 1)
                    name = 'place-' + str(mol_idx + 1) + '-' + str(i + 1)
                    mol, ddg, rmsd = place_smiles(frag_mols, smi, protein, work_dir=work_dir, name=name, seed=None)
                    if not mol:
                        DmLog.emit_event("Placement failed for molecule", count, smiles)
                    else:
                        num_placements += 1
                        mol = prepare_mol_for_writing(mol, num_placements, id, protein, frag_mols, ddg, rmsd, removeHydrogens=removeHydrogens,
                                                      idxProp=idxProp, ddgProp=ddgProp, rmsdProp=rmsdProp, smilesProp=smilesProp,
                                                      fragIdField=fragIdField, refMolsProp=refMolsProp,
                                                      proteinProp=proteinProp, proteinPropValue=proteinPropValue)
                        writer.write(mol)
                mol_idx += 1

            except Exception as e:
                utils.log("Failed to place molecule", count, e)
                print(traceback.format_exc())
                num_errors += 1

        if count == num_errors: # all attempts fail
            DmLog.emit_event('All invocations failed')
            raise ValueError('All invocations failed')

    return len(smiles), num_placements, num_errors


def combine_fragments(fragments, protein, outputFileName, work_dir=None, fragIdField=None, removeHydrogens=True,
                      refMolsProp=None, idxProp='IDX', ddgProp='DDG', rmsdProp='RMSD', smilesProp=None,
                      proteinProp=None, proteinPropValue=None, num_runs=1):
    """
    Combine these fragments into a single monster molecule

    :param fragments:
    :param protein:
    :param outputFileName:
    :param work_dir:
    :param fragIdField:
    :param refMolsProp:
    :param removeHydrogens:
    :param idxProp:
    :param ddgProp:
    :param rmsdProp:
    :param smilesProp:
    :param proteinProp:
    :param proteinPropValue:
    :param num_runs:
    :return:
    """

    if not os.path.isfile(protein):
        raise ValueError("Protein file not found: " + protein)

    num_merges = 0
    num_failures = 0
    num_errors = 0

    frag_mols = utils.read_molecules(fragments)
    DmLog.emit_event('Read', len(frag_mols), 'fragments')

    utils.expand_path(outputFileName)

    with Chem.SDWriter(outputFileName) as writer:
        for i in range(num_runs):
            DmLog.emit_event('Running combination', i + 1)
            try:
                v = create_victor(frag_mols, protein, work_dir, None)
                v.combine(long_name='combine-' + str(i + 1))
                mol = v.minimized_mol
                if not mol:
                    num_failures += 1
                else:
                    num_merges += 1
                    mol = prepare_mol_for_writing(mol, num_merges, num_merges, protein, frag_mols, v.ddG, v.mrmsd.mrmsd,
                                                  idxProp=idxProp, ddgProp=ddgProp, rmsdProp=rmsdProp,
                                                  smilesProp=smilesProp, removeHydrogens=removeHydrogens,
                                                  refMolsProp=refMolsProp, fragIdField=fragIdField,
                                                  proteinProp=proteinProp, proteinPropValue=proteinPropValue)
                    writer.write(mol)

            except Exception as e:
                DmLog.emit_event("Failed to combine molecules", i, e)
                utils.log(traceback.format_exc())
                num_errors += 1

    if num_runs == num_errors: # all attempts fail
        DmLog.emit_event('All invocations failed')
        raise ValueError('All invocations failed')

    return num_merges, num_failures, num_errors


def main():

    # Example (combine):
    #   ./merger.py -f data/diamond-x0216_A.mol data/diamond-x0259_A.mol -p data/mac1.pdb -o output/foo.sdf
    #
    # Example (placement):
    #   ./merger.py -f data/x0032_0A.mol data/x0103_0A.mol -p data/apo_example1.pdb --smiles-file data/inputs.smi -o output/foo.sdf
    #
    #   ./merger.py -f data/x0032_0A.mol data/x0103_0A.mol -p data/apo_example1.pdb --smiles-strings 'COC(=O)c1[nH]c(C)c(C(=O)C)c1C' 'C=CCSc1nnc(NCc2ccccc2)s1' -o output/foo.sdf

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Fragmenstein place and combine')
    parser.add_argument('-f', '--fragments', nargs='+', required=True, help="Molfiles or SD-files with fragments")
    parser.add_argument('--frag-id-field', help="Field with the fragment IDs")
    parser.add_argument('--smiles-strings', nargs="+", help="SMILES to place (if not supplied then fragments are combined)")
    parser.add_argument('--smiles-file',
                        help="SMILES file to place (if not supplied then fragments are combined). SMILES is assumed to be in first column")
    parser.add_argument('-d', '--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--smiles-id-field', type=int, default=1, help="Column for the SMILES ID field (zero based integer)")
    parser.add_argument('--read-header', action='store_true',
                        help="Read a header line with the field names when reading smiles")
    parser.add_argument('-p', '--protein', required=True, help="PDB file for protein")
    parser.add_argument('-o', '--outfile', required=True, help="Output file (.sdf)")
    parser.add_argument('-w', '--work-dir', help="Location of the work dir. If not defined then ./output is used")
    parser.add_argument('-k', '--keep-hydrogens', action='store_true', help='Keep hydrogens in the outputs')
    parser.add_argument('-c', '--count', type=int, default=1,
                        help="Number of structures to generate for each set of inputs")
    parser.add_argument('--idx-prop-name', help='Name for the IDX property (default IDX)')
    parser.add_argument('--ddg-prop-name', help='Name for the ddG property (default DDG)')
    parser.add_argument('--rmsd-prop-name', help='Name for the RMSD property (default RMSD)')
    parser.add_argument('--ref-mols-prop-name', help='Name for the field with the input IDs')
    parser.add_argument('--smiles-prop-name', help='Name for the field containing the SMILES of the generated molecule')
    parser.add_argument('--protein-prop-name', help='Write the protein file name as this property')
    parser.add_argument('--protein-prop-value', help='Value to use for the protein name. If not specified the file name is used')

    args = parser.parse_args()
    DmLog.emit_event("merger: ", args)

    if args.smiles_file and args.smiles_strings:
        raise ValueError("Can't specify both smiles-file and smiles-strings arguments.")

    delimiter = utils.read_delimiter(args.delimiter)

    kwargs = {'num_runs': args.count}
    if args.work_dir:
        kwargs['work_dir'] = args.work_dir
    if args.keep_hydrogens:
        kwargs['removeHydrogens'] = not args.keep_hydrogens
    if args.idx_prop_name:
        kwargs['idxProp'] = args.idx_prop_name
    if args.ddg_prop_name:
        kwargs['ddgProp'] = args.ddg_prop_name
    if args.rmsd_prop_name:
        kwargs['rmsdProp'] = args.rmsd_prop_name
    if args.ref_mols_prop_name:
        kwargs['refMolsProp'] = args.ref_mols_prop_name
    if args.frag_id_field:
        kwargs['fragIdField'] = args.frag_id_field
    if args.smiles_prop_name:
        kwargs['smilesProp'] = args.smiles_prop_name
    if args.protein_prop_name:
        kwargs['proteinProp'] = args.protein_prop_name
    if args.protein_prop_value:
        kwargs['proteinPropValue'] = args.protein_prop_value

    t0 = time.time()
    if args.smiles_file:
        smiles, ids = read_smiles_file(args.smiles_file, delimiter=delimiter, readHeader=args.read_header, idColumn=args.smiles_id_field)
        num_mols, num_placements, num_failures = \
            place_multiple_smiles(args.fragments, smiles, ids, args.protein, args.outfile, **kwargs)

        t1 = time.time()
        # Duration? No less than 1 second?
        duration_s = int(t1 - t0)

        DmLog.emit_event('Placing {} smiles generated {} structures, {} failures. Took {}s'.format(
            num_mols, num_placements, num_failures, duration_s))
        DmLog.emit_cost(num_placements)

    elif args.smiles_strings:
        ids = []
        smiles = []
        for i, smi in enumerate(args.smiles_strings):
            smiles.append(smi)
            ids.append(str(i + 1))

        num_mols, num_placements, num_failures = \
            place_multiple_smiles(args.fragments, smiles, ids, args.protein, args.outfile, **kwargs)

        t1 = time.time()
        # Duration? No less than 1 second?
        duration_s = int(t1 - t0)

        DmLog.emit_event('Placing {} smiles generated {} structures, {} failures. Took {}s'.format(
            num_mols, num_placements, num_failures, duration_s))
        DmLog.emit_cost(num_placements)

    else:
        num_merges, num_failures, num_errors = combine_fragments(args.fragments, args.protein, args.outfile, **kwargs)
        t1 = time.time()
        # Duration? No less than 1 second?
        duration_s = int(t1 - t0)

        DmLog.emit_event('Combining {} molecules generated {} structures, {} failures, {} errors. Took {}s'.format(
            len(args.fragments), num_merges, num_failures, num_errors, duration_s))
        DmLog.emit_cost(num_merges)


if __name__ == "__main__":
    main()
