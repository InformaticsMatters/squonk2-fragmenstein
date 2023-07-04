This describes how to run the `fragmenstein-combine`, `fragmenstein-place-file` and `fragmenstein-place-string` jobs 
from the `fragmenstein` category in the `fragmenstein` collection.

## What the job does

This job uses the fragmenstein tools to generate molecules that are a merge of 2 fragments (combine) or place molecules 
(as SMILES) onto the 3D coordinates of fragments. These "monster" molecules contain features from the input molecules
and aim to represent a merge of the features of the input molecules. Generated molecules are minimised within the
active site of the corresponding protein.

It may not be possible to generate molecules for all combinations of fragments and some molecules that are generated
will not be chemically feasible. See the
[Fragmenstein documentation](https://github.com/matteoferla/Fragmenstein) for more information about Fragmenstein.

## Implementation details

* Job implementation: [merger.py](/merger.py)
* Job definition: `jobs.fragmenstein-combine` in [fragmenstein.yaml](../fragmenstein.yaml) (combine fragments)
* Job definition: `jobs.fragmenstein-place-file` in [fragmenstein.yaml](../fragmenstein.yaml) (place molecule onto fragments)
* Job definition: `jobs.fragmenstein-place-string` in [fragmenstein.yaml](../fragmenstein.yaml) (place molecule onto fragments)


## How to run the job

### Inputs

* **Fragment molecules**: Two or more fragment molecules as molfile or SD-file (1).
* **PDB file for protein**: The protein in which the generated molecules are minimised.
* **File with SMILES to place**: File with smiles to place (SMILES must be the first column). Only needed for the
  `jobs.fragmenstein-place-file` mode.

Notes:
(1) Input fragments are specified either as separate molfiles or a single SD-file. For the place and combine modes
specify 2 or 3 input fragments.

### Options

All modes:

* **Output file name**: Name for the output SD-file.
* **Number of molecules to generate**: Number of molecules to generate for each combination of fragments/inputs.
* **Keep hydrogens in the outputs**: Retain hydrogens in the output. By default they are removed.
* **Input field name containing the fragment ID**: Field in the input fragments that will be written to the outputs as
the `ref_mols` field.
* **Include SMILES in output using this field name**: If specified the SMILES of the generated molecule is written to
a field of this name in the output molecules.
* **Include PDB details in output using this field name**: If specified details of the PDB file are written to a field
of this name in the output molecules.
* **Use this value for the proteinFieldName**: Use this value for the proteinFieldName value. If not specified the PDB file name is used.

### Outputs

A SD-file is output containing 3D structures for the generated molecules. *Number of molecules to generate* molecules are
generated for each set of input fragments. Some input fragments may not generate any outputs.

The following fields are written to the SD-file:

* **IDX**: The index of the molecule (related to *Number of molecules to generate*).
* **DDG**: The difference in delta G between the bound and unbound minimised molecule.
* **RMSD** The mean RMSD between the fragments and the minimised molecule.
* **Value of "Includes SMILES in output using this field name"**: If this option is specified the SMILES of the generated
  molecule is written to this field
**Value of "Include PDB details in output using this field name"**: If this option is specified details of the PDB file 
are written to this field. The value written is the value of the "Use this value for the proteinFieldName" property if
defined, or the PDB file name if not.
* **ref_mols**: If the *Input field name containing the fragment ID* option is specified then the value of this field is
used as the ID of the input fragment and these IDs are written to this field as a comma separated list (2).

Notes:
(2) e.g. if the input fragments have their IDs as their title line (first line in the record) then you can specify 
`_Name` as the value for this field and those IDs will be written to the `ref_mols` field e.g. as `Mpro-x0072_0A,Mpro-x0104_0A`.
If the input is a SD-file then you can also specify the name of any regular SD-file field. If you do not specify a value
for this option then the `ref_mols` field is not written.
