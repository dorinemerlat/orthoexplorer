#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ORTHOEXPLORER } from './workflows/orthoexplorer'

workflow {
    ORTHOEXPLORER()
}