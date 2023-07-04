import os
import sys
from rdkit import Chem


def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)


def read_delimiter(input):
    if input:
        if 'tab' == input:
            delimiter = '\t'
        elif 'space' == input:
            delimiter = None
        elif 'comma' == input:
            delimiter = ','
        elif 'pipe' == input:
            delimiter = '|'
        else:
            delimiter = input
    else:
        delimiter = None
    return delimiter


def expand_path(path):
    """
    Create any necessary directories to ensure that the file path is valid

    :param path: a filename or directory that might or not exist
    """
    head_tail = os.path.split(path)
    if head_tail[0]:
        if not os.path.isdir(head_tail[0]):
            log('Creating directories for', head_tail[0])
            os.makedirs(head_tail[0], exist_ok=True)

def read_molecules(inputs):
    """
    Read input molecules.
    A list of inputs is specified, each one containing either a single filename or a comma separated list of filenames
    (no spaces). The files can either be .mol files with a single molecule or .sdf files with multiple molecules.
    The molecules in all the files specified are read and returned as a list of molecules.
    Examples:
    ['mols.sdf']                          - single SDF with one or more molecules
    ['mol1.mol', 'mol2.mol', 'mol3.mol']  - 3 molfiles as separate elements of the list
    ['mol1.mol,mol2.mol,mol3.mol']        - 3 molfiles as a single comma separated element of the list
    ['mol1.mol,mol2.mol', 'mols.sdf']     - 3 molfiles plus 1 SDF

    :param inputs: Input filenames
    :return: List of molecules that have been read
    """
    mols = []
    for input in inputs:
        tokens = input.split(',')
        for token in tokens:
            if token.endswith('.mol'):
                m = Chem.MolFromMolFile(token)
                if m:
                    mols.append(m)
                else:
                    DmLog.emit_event('WARNING: could not process', token)
            elif token.endswith('.sdf'):
                supplr = Chem.ForwardSDMolSupplier(token)
                idx = 0
                for m in supplr:
                    if m:
                        mols.append(m)
                    else:
                        DmLog.emit_event('WARNING: could not process molecule', idx, 'from', token)
            else:
                raise ValueError("Invalid file", token)
    return mols