# OrthoExplorer

OrthoExplorer is a modular Nextflow DSL2 pipeline for preparing proteomes, inferring orthogroups with OrthoFinder, and performing downstream gene-conservation and pangenome analyses.

## Workflow

The pipeline is divided into three subworkflows:

1. **Dataset preparation**
   - validates species and outgroup tables;
   - uses user-provided proteomes when available;
   - otherwise downloads a RefSeq genome and GFF3 from NCBI;
   - extracts the longest protein isoform for each gene with AGAT;
   - downloads the NCBI taxonomy for all taxa.

2. **Orthology inference**
   - collects all prepared proteomes;
   - runs OrthoFinder;
   - exports normalized OrthoFinder result filenames.

3. **Downstream analyses**
   - reformats OrthoFinder species names;
   - assigns gene-conservation ranks;
   - defines pangenome categories;
   - generates Venn/UpSet summaries.

## Input tables

Both `species.csv` and `outgroups.csv` must contain the following columns:

```csv
taxid,name,file
126957,Strigamia maritima,/path/to/strigamia-maritima.faa
1490507,Geophilus flavus,
```

`taxid` and `name` are mandatory and must be non-empty.

`file` is optional:

- when populated, the referenced protein FASTA is used;
- when empty, OrthoExplorer downloads a RefSeq genome and GFF3 and extracts one canonical protein isoform per gene.

Paths in `file` should be absolute, or valid relative to the directory from which Nextflow is launched.

## Required software

The current modules expect these commands to be available through the execution environment:

- Nextflow;
- NCBI Datasets CLI;
- `unzip`;
- AGAT;
- OrthoFinder;
- `download_taxonomy.py`;
- `reformat_orthofinder_tables.py` (included in `bin/`);
- `assign_gene_conservation_rank.py`;
- `venn_upset_orthogroups.py`.

Software environments should preferably be defined with containers or Conda in the Nextflow configuration rather than with hard-coded `module load` commands inside process modules.

## Configuration

Example configuration:

```groovy
params {
    species   = "${projectDir}/species_myriapoda.csv"
    outgroups = "${projectDir}/outgroups_myriapoda.csv"
    groups    = "Diplopoda,Chilopoda,Hexapoda,Crustacea,Chelicerata"
    colors    = "${projectDir}/data/colors.yaml"

    outdir    = "${projectDir}/results"
    cache_dir = "${projectDir}/cache/work"
    conda_env = "${projectDir}/envs/orthoexplorer"
}

workDir = params.cache_dir
```

Process resources are defined in `conf/modules.config`, while publication rules are defined in `conf/publish.config`.

## Running the pipeline

```bash
nextflow run main.nf \
    -c nextflow_myriapoda.config \
    -resume
```

Parameters may also be overridden from the command line:

```bash
nextflow run main.nf \
    -c nextflow_myriapoda.config \
    --species species_myriapoda.csv \
    --outgroups outgroups_myriapoda.csv \
    --colors data/colors.yaml \
    --outdir results \
    -resume
```

## Stub testing

Each process contains a `stub` block that checks whether its main command is available and creates the expected empty outputs.

```bash
nextflow run main.nf \
    -c nextflow_myriapoda.config \
    -stub-run
```

The stub run validates workflow wiring and output declarations, but it does not validate the biological content of the outputs.

## Output structure

```text
results/
├── taxonomy/
├── orthofinder/
├── gene_conservation/
└── pangenome/
```

The Nextflow execution cache is stored separately under `cache/work/`.

---


## Credits

Originally written by Dorine Merlat (dorine.merlat@etu.unistra.fr).
Thanks to Arnaud Kress and Odile Lecompte for assistance.


# Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

---

# License

This project is distributed under the MIT License.

