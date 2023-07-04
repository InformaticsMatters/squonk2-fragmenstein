params.inputs = 'inputs.sdf'

params.min_dist_thr = 3
params.min_num = 2
params.max_num = 2 // or 3
params.min_dist = 1.5


process pairwise_prep {

    container 'informaticsmatters/squonk2-fragmenstein:stable'

    input:
    path inputs

    output:
    path '*.sdf'

    """
    /code/prep_compatible_frags.py -i '${inputs.join("\' \'")}' --min-num $params.min_num --max-num $params.max_num --min-dist $params.min_dist
    """
}