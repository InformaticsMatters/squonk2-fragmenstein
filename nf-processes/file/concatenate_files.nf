params.publish_dir = ''
params.publish_dir_mode = 'copy'
params.optional = false

/** Concatenate the files matching glob into outputfile.
*/
process concatenate_files {

    container 'informaticsmatters/vs-prep:3.1.0'
    // 'publishDir' inside an 'if' is not accepted by the Nextflow 25+ parser;
    // 'enabled:' is the supported form. The elvis guards the path, which must
    // still be a valid string even when publishing is disabled.
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path part
    val outputfile // e.g. 'results.sdf'
    val glob       // e.g. 'merged_*.sdf'

    output:
    path "${outputfile}", optional: params.optional

    script:
    """
    DIR=\$(dirname "${outputfile}")
    mkdir -p \$DIR
    ls ${glob} | xargs cat >> ${outputfile}
    """
}
