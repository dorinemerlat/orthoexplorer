#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ORTHOEXPLORER } from './workflows/orthoexplorer'


def format_dataset_row(row, dataset_type, source_csv) {

    def taxid = row.taxid?.toString()?.trim()
    def name  = row.name?.toString()?.trim()
    def fasta = row.fasta?.toString()?.trim()
    def gff   = row.gff?.toString()?.trim()

    if (!taxid || !name) {
        error """
        Invalid row in '${source_csv}'.

        Columns 'taxid' and 'name' are mandatory.

        Offending row:
        ${row}
        """.stripIndent()
    }

    /*
     * Protein FASTA and GFF must either both be supplied
     * or both be left empty.
     */
    if ((fasta && !gff) || (!fasta && gff)) {
        error """
        Invalid row for '${name}' in '${source_csv}'.

        Columns 'fasta' and 'gff' must either both contain a file
        or both be empty.

        fasta = ${fasta ?: '<empty>'}
        gff   = ${gff   ?: '<empty>'}
        """.stripIndent()
    }

    def id = name
        .toLowerCase()
        .replaceAll(/[^a-z0-9]/, '-')
        .replaceAll(/^-|-$/, '')

    def meta = [
        id           : id,
        name         : name,
        taxid        : taxid,
        dataset_type : dataset_type
    ]

    return [id, meta, fasta ?: null, gff   ?: null]
}


workflow {

    if (!params.ingroups) {
        error "Missing required parameter: --ingroups"
    }

    if (!params.outgroups) {
        error "Missing required parameter: --outgroups"
    }


    ingroups_csv = file(params.ingroups, checkIfExists: true)
    outgroups_csv = file(params.outgroups, checkIfExists: true)
    tree = params.tree ? file(params.tree, checkIfExists: true) : null
    colors_yaml = params.colors ? file(params.colors, checkIfExists: true) : null

    /*
     * Read and validate the two input CSV files.
     */
    ingroups = Channel
        .fromPath(ingroups_csv, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            format_dataset_row(row, 'ingroup', ingroups_csv)
        }

    outgroups = Channel
        .fromPath(outgroups_csv, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            format_dataset_row(row, 'outgroup', outgroups_csv)
        }
    
    tree = Channel.fromPath(params.tree, checkIfExists: true)

    colors = params.colors ? Channel.fromPath(params.colors, checkIfExists: true) : null

    all_datasets = ingroups
        .concat(outgroups)

    /*
     * Separate datasets with supplied annotations from datasets
     * that must be downloaded.
     */
    all_datasets
        .branch {
            supplied:
                it[2] && it[3]

            download:
                !it[2] && !it[3]
        }
        .set { input_datasets }


    /*
     * Print summary once before the workflow is launched.
     */
    ingroups_count = ingroups.count()
    outgroups_count = outgroups.count()
    all_datasets_count = all_datasets.count()
    supplied_count = input_datasets.supplied.count()
    download_count = input_datasets.download.count()

    summary_counts = ingroups_count
        .combine(outgroups_count)
        .combine(supplied_count)
        .combine(download_count)
        .combine(all_datasets_count)

    summary_counts.view { ingroup_n, outgroup_n, supplied_n, download_n, total_n ->

        """
        ============================================================
        OrthoExplorer
        ============================================================
        Ingroups table       : ${ingroups_csv}
        Outgroups table      : ${outgroups_csv}
        Tree                 : ${params.tree ?: 'OrthoFinder species tree'}
        Groups               : ${params.clades ?: 'not specified'}
        Colors               : ${colors_yaml ?: 'not specified'}

        Ingroups             : ${ingroup_n}
        Outgroups            : ${outgroup_n}
        Supplied annotations : ${supplied_n}
        To download          : ${download_n}
        Total datasets       : ${total_n}

        Output directory     : ${params.outdir}
        Work directory       : ${workflow.workDir}
        Launch directory     : ${workflow.launchDir}
        Profile              : ${workflow.profile ?: 'standard'}
        ============================================================
        """.stripIndent()
    }

    /*
     * Convert user-provided paths into staged Nextflow files.
     */
    input_datasets.supplied
        .map { id, meta, fasta, gff ->
            [ id, meta, file(fasta, checkIfExists: true), file(gff, checkIfExists: true) ]
        }
        .set { supplied_datasets }

    /*
     * Only metadata are required by the download process.
     */
    input_datasets.download
        .map { id, meta, fasta, gff ->
            [ id, meta ]
        }
        .set { datasets_to_download }
    
    ORTHOEXPLORER(
        supplied_datasets,
        datasets_to_download,
        ingroups_csv,
        outgroups_csv,
        tree,
        colors
    )
}