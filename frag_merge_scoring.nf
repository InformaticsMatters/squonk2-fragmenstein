/* Copyright 2022 Informatics Matters Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*
Run fragmenstein combine of 2-wise or 3-wise combinations of fragments

Example:
nextflow run frag_merge_scoring.nf --fragments 'data/Mpro-x0072_0A.mol,data/Mpro-x0104_0A.mol,data/Mpro-x0107_0A.mol' \
  --protein 'data/Mpro-x1249_0A_apo-desolv.pdb' \
  --outfile merged.sdf
*/

nextflow.enable.dsl=2

params.scratch = false
params.fragments = 'fragments.sdf'
params.protein = 'protein.pdb'
params.outfile = 'merges.sdf'
params.publish_dir = './'

// files - specified as comma separated list of files as a single string (no spaces) to the --fragments argument
fragments = Channel.of(params.fragments.toString())
                    .splitCsv()
                    .flatten()
                    .map {it -> file(it, checkIfExists:true)}
                    .toList()
protein = file(params.protein)

// includes
include { pairwise_prep } from './nf-processes/fragmenstein/prep_compatible_frags.nf'
include { combine } from './nf-processes/fragmenstein/fragmenstein_combine.nf'
include { scoring } from './nf-processes/xchem/scoring.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf' addParams(
    outputfile: params.outfile,
    glob: 'scored_*.sdf')

dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))
def curr_t() { dateFormat.format(new java.util.Date()) }

def wrkflw = 'fragmenstein_combine_scoring'
now = curr_t()
log.info("$now # PROGRESS -START- $wrkflw:pairwise_prep 1")

// workflows
workflow fragmenstein_combine_scoring {

    take:
    fragments
    protein

    main:
    pairwise_prep(fragments)
    combine(pairwise_prep.out.flatten(), protein)
    scoring(combine.out[0], combine.out[1])
    concatenate_files(scoring.out[0].collect())

    int cost = 0
    int combination_count = 0
    int fragmenstein_count = 0
    int scoring_count = 0

    pairwise_prep.out.flatten().subscribe {
        now = curr_t()
        if (combination_count == 0) log.info("$now # PROGRESS -DONE- $wrkflw:pairwise_prep 1")
        combination_count++
        log.info("$now # PROGRESS -START- $wrkflw:combine $combination_count")
    }

    combine.out[2].subscribe {
        cost += new Integer(it)
        fragmenstein_count += 1
        now = curr_t()
        log.info("$now # INFO -COST- $cost $fragmenstein_count")
        log.info("$now # PROGRESS -DONE- $wrkflw:combine $fragmenstein_count")
        log.info("$now # PROGRESS -START- $wrkflw:scoring $fragmenstein_count")
    }

    scoring.out[1].subscribe {
        cost = new Integer(it)
        scoring_count++
        now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:scoring $scoring_count")
    }

    scoring.out[0].collect().subscribe {
        now = curr_t()
        log.info("$now # PROGRESS -START- $wrkflw:concatenate_files 1")
    }

    concatenate_files.out.subscribe {
        now = curr_t()
        log.info("$now # PROGRESS -DONE- $wrkflw:concatenate_files 1")
    }

    emit:
    concatenate_files.out
}

workflow {
    fragmenstein_combine_scoring(fragments, protein)
}
