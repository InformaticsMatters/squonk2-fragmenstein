This describes how to run the `fragmenstein-combine-multi` and `fragmenstein-combine-multi-scoring` jobs from the `fragmenstein`
category in the `fragmenstein` collection.

## What the job does

Like the simple `fragmenstein-combine` job, this job uses the fragmenstein tools to generate molecules that are a merge
of 2 or more fragments. This variant is a Nextflow workflow that takes multiple input fragments, generates all possible
pairwise (or 3-wise) combinations of these and runs fragmenstein on those different subsets of fragments.

A naive approach might be to use all the fragments obtained in a fragment screen as inputs.

Following the multiple independent executions of fragmenstein the resulting molecules are concatenated into a single results
file. For the `fragmenstein-combine-multi-scoring` variant these additional scoring algorithms are applied to the resulting
molecules to assist in assessing their goodness:

* SA Score - a synthetic accessibility score
* SuCOS - a score of feature and volume overlap with the source fragments
* Open3DAlign - a score of the rigid 3D alignment of molecule onto the source fragments

See the [Fragmenstein documentation](https://github.com/matteoferla/Fragmenstein) for more information about Fragmenstein.


## Implementation details

* Job implementation: [fragmenstein-combine-multi](/frag_merge.nf), [fragmenstein-combine-multi-scoring](/frag_merge_scoring.nf)
* Job definition: `jobs.fragmenstein-combine-multi` in [fragmenstein.yaml](../fragmenstein.yaml),
  `jobs.fragmenstein-combine-multi-scoring` in [fragmenstein.yaml](../fragmenstein.yaml)

For details about the individual steps see these sources:

* Generation of fragment combinations: [Python module](/prep_compatible_frags.py)
* Fragmenstein [Job docs](fragmenstein.md), [Python module](/merger.py)
* SA Score [Job docs](https://github.com/InformaticsMatters/virtual-screening/blob/main/data-manager/docs/rdkit/sa-score.md),
  [Python module](https://github.com/InformaticsMatters/virtual-screening/blob/main/sa_score.py)
* SuCOS score [Job docs](https://github.com/InformaticsMatters/virtual-screening/blob/main/data-manager/docs/xchem/sucos.md),
  [Python module](https://github.com/InformaticsMatters/virtual-screening/blob/main/sucos.py)
* Open3DAlign score [Job docs](https://github.com/InformaticsMatters/virtual-screening/blob/main/data-manager/docs/rdkit/open3dalign.md),
  [Python module](https://github.com/InformaticsMatters/virtual-screening/blob/main/open3dalign.py)

## How to run the job

### Inputs

* **Fragment molecules**: Two or more fragment molecules as molfile or SD-file (1).
* **PDB file for protein**: The protein in which the generated molecules are minimised.

Notes:
(1) Input fragments are specified either as separate molfiles or a single SD-file. Specify a relatively small number of
fragments (e.g. less than 20) and compatible 2-way and 3-way combinations of these are evaluated.

### Options

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
* **Minimum number to combine**: The minimum number of fragments to try to combine (2 or 3) (2).
* **Maximum number to combine**: The maximum number of fragments to try to combine (2 or 3) (2).
* **Maximum distance**: Do not try to combine fragments that are greater than this distance apart (default 1.5A).

Notes:
(2) Typically use 2 for min and max, but also consider using 3 for max.

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
used as the ID of the input fragment and these IDs are written to this field as a comma separated list (3).
* **sa_score**: The synthetic accessibility score. 1 means easy, 10 means hard.
* **SuCOS_Score**: SuCOS score (0=none, 1=perfect)
* **SuCOS_FeatureMap_Score**: SuCOS feature map score (0=none, 1=perfect)
* **SuCOS_Protrude_Score**: SuCOS protrude score (0=none, 1=perfect)
* **o3da_score**: Open3DAlign score (larger numbers indicate better alignment)
* **o3da_score_rel**: Open3DAlign relative score (0=none, 1=perfect)
* **o3da_align**: Open3DAlign alignment score (larger numbers indicate better alignment)

Notes:
(3) e.g. if the input fragments have their IDs as their title line (first line in the record) then you can specify 
`_Name` as the value for this field and those IDs will be written to the `ref_mols` field e.g. as `Mpro-x0072_0A,Mpro-x0104_0A`.
If the input is a SD-file then you can also specify the name of any regular SD-file field. If you do not specify a value
for this option then the `ref_mols` field is not written.
