params.crippen = false
params.remove_hydrogens = true
params.score_mode = 'all'
params.tanimoto = false

process scoring {

    container 'informaticsmatters/vs-prep:latest'

    input:
    path inputs // .sdf
    path query  // .sdf or .mol

    output:
    path "scored_${inputs.name}"
    env COUNT

    """
    OUT='scored_${inputs.name}'

    /code/sa_score.py -i '$inputs' -o tmp1.sdf

    /code/sucos.py\
      --input tmp1.sdf\
      --reference $query\
      --output tmp2.sdf\
      --score-mode '$params.score_mode'\
      ${params.tanimoto ? '--tanimoto' : ''}

    /code/open3dalign.py\
      --inputs tmp2.sdf\
      --query $query\
      --outfile "\$OUT"\
      ${params.remove_hydrogens ? '--remove-hydrogens' : ''}\
      ${params.crippen ? '--crippen' : ''}

    # count the number of outputs - for some strange reason the fgrep command fails is the file is empty
    if [ -s "\$OUT" ]
    then
      COUNT=\$(fgrep -c '\$\$\$\$' "\$OUT")
    else
      COUNT=0
    fi
    """
}
