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

// includes
include { pairwise_prep } from './nf-processes/fragmenstein/prep_compatible_frags.nf'
include { combine } from './nf-processes/fragmenstein/fragmenstein_combine.nf'
include { scoring } from './nf-processes/xchem/scoring.nf'
include { concatenate_files } from './nf-processes/file/concatenate_files.nf'

// Self-contained: the parser does not allow statements at the top level, so
// the formatter is built per call rather than held in a script-level variable.
def curr_t() {
    def dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+00:00'", Locale.UK)
    dateFormat.setTimeZone(TimeZone.getTimeZone('UTC'))
    return dateFormat.format(new java.util.Date())
}

// workflows
workflow fragmenstein_combine_scoring {

    take:
    fragments
    protein

    main:
    def wrkflw = 'fragmenstein_combine_scoring'
    log.info("${curr_t()} # PROGRESS -START- $wrkflw:pairwise_prep 1")

    pairwise_prep(fragments)
    combine(pairwise_prep.out.flatten(), protein)
    scoring(combine.out[0], combine.out[1])
    concatenate_files(scoring.out[0].collect(), params.outfile, 'scored_*.sdf')

    // Counters are AtomicInteger rather than int: the parser rejects '++',
    // and these subscribe callbacks can run on different threads.
    def cost = new java.util.concurrent.atomic.AtomicInteger()
    def combination_count = new java.util.concurrent.atomic.AtomicInteger()
    def fragmenstein_count = new java.util.concurrent.atomic.AtomicInteger()
    def scoring_count = new java.util.concurrent.atomic.AtomicInteger()

    pairwise_prep.out.flatten().subscribe { _part ->
        def now = curr_t()
        if (combination_count.get() == 0) log.info("$now # PROGRESS -DONE- $wrkflw:pairwise_prep 1")
        log.info("$now # PROGRESS -START- $wrkflw:combine ${combination_count.incrementAndGet()}")
    }

    combine.out[2].subscribe { count_file ->
        def total = cost.addAndGet(count_file.text.trim() as Integer)
        def n = fragmenstein_count.incrementAndGet()
        def now = curr_t()
        log.info("$now # INFO -COST- $total $n")
        log.info("$now # PROGRESS -DONE- $wrkflw:combine $n")
        log.info("$now # PROGRESS -START- $wrkflw:scoring $n")
    }

    // NOTE: this deliberately preserves the pre-existing behaviour, which
    // *sets* the running cost rather than adding to it, and emits no -COST-
    // line of its own. See squonk2-jobs#72 - changing it changes what is
    // billed, so it is not being changed here.
    scoring.out[1].subscribe { count_file ->
        cost.set(count_file.text.trim() as Integer)
        def n = scoring_count.incrementAndGet()
        log.info("${curr_t()} # PROGRESS -DONE- $wrkflw:scoring $n")
    }

    scoring.out[0].collect().subscribe { _results ->
        log.info("${curr_t()} # PROGRESS -START- $wrkflw:concatenate_files 1")
    }

    concatenate_files.out.subscribe { _result ->
        log.info("${curr_t()} # PROGRESS -DONE- $wrkflw:concatenate_files 1")
    }

    emit:
    concatenate_files.out
}

workflow {
    // files - specified as comma separated list of files as a single string
    // (no spaces) to the --fragments argument
    def fragments = Channel.of(params.fragments.toString())
                        .splitCsv()
                        .flatten()
                        .map { it -> file(it, checkIfExists: true) }
                        .toList()
    def protein = file(params.protein)

    fragmenstein_combine_scoring(fragments, protein)
}
