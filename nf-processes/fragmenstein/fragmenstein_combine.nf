params.inputs = 'inputs.sdf'
params.protein = 'protein.pdb'
params.count = 1
params.keephs = false
params.ref_mols_prop = null
params.frag_id_field = null
params.smiles_prop = null
params.protein_prop_name = null
params.protein_prop_value = null


process combine {

    container 'informaticsmatters/squonk2-fragmenstein:stable'

    input:
    path inputs
    path protein

    output:
    path 'merged_*.sdf'
    path 'fragments.sdf'
    env COUNT

    """
    /code/merger.py -f '$inputs' -p '$protein' -o 'merged_$inputs' ${params.keephs ? '--keep-hydrogens' : ''}\
      --count $params.count\
      ${params.ref_mols_prop ? '--ref-mols-prop-name \'' + params.ref_mols_prop + '\'' : ''}\
      ${params.frag_id_field ? '--frag-id-field \'' + params.frag_id_field + '\'' : ''}\
      ${params.smiles_prop ? '--smiles-prop-name \'' + params.smiles_prop + '\'' : ''}\
      ${params.protein_prop_name ? '--protein-prop-name \'' + params.protein_prop_name + '\'' : ''}\
      ${params.protein_prop_value ? '--protein-prop-value \'' + params.protein_prop_value + '\'' : ''}

      # downstream processes need to use the fragments
      cp '$inputs' fragments.sdf

      # record the number of outputs
      COUNT=$params.count
    """
}