process MOSDEPTH_REGIONS_TO_MQC {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.8.3'
        : 'biocontainers/python:3.8.3'}"

    input:
    tuple val(meta), path(regions_bed), path(thresholds_bed)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def header = [
        "# id: '${prefix}_mosdepth_per_region'",
        "# section_name: 'Per-region coverage'",
        "# description: 'Mean coverage and bases covered at ≥20X and ≥30X per target region'",
        "# order: 600",
        "# plot_type: 'table'",
        "# pconfig:",
        "#    id: '${prefix}_mosdepth_per_region_table'",
        "Gene\\tChrom\\tStart\\tEnd\\tMean coverage\\t% bases ≥20X\\t% bases ≥30X",
    ].join("\\n")
    """
    echo -e "${header}" > ${prefix}_mqc.tsv
    python3 -c "
import gzip
def read_bed(path):
    with gzip.open(path, 'rt') as f:
        return [line.strip().split('\\t') for line in f if line.strip() and not line.startswith('#')]
regions = read_bed('${regions_bed}')
thresholds = read_bed('${thresholds_bed}')
for r, t in zip(regions, thresholds):
    length = int(r[2]) - int(r[1])
    pct_20x = round(int(t[4]) / length * 100, 2) if length > 0 else 0
    pct_30x = round(int(t[5]) / length * 100, 2) if length > 0 else 0
    print(r[3], r[0], r[1], r[2], r[4], pct_20x, pct_30x, sep='\\t')
" >> ${prefix}_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mqc.tsv
    """
}
