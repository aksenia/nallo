process CRAMINO_TO_MQC {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(cramino_txt)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def section_name = task.ext.section_name ?: 'Read-level QC'
    def header = [
        "# id: '${prefix}_cramino'",
        "# section_name: '${section_name}'",
        "# order: 500",
        "# plot_type: 'table'",
        "# pconfig:",
        "#    id: '${prefix}_cramino_table'",
        "Sample\\tNumber of reads\\tYield, Gb\\tYield, Gb (>25kb)\\tN50\\tMedian length\\tMean length\\tMedian identity",
    ].join("\\n")
    """
    echo -e "${header}" > ${prefix}_mqc.tsv
    python3 -c "
import sys
stats = {}
with open('${cramino_txt}') as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) == 2:
            stats[parts[0]] = parts[1]
row = [
    '${prefix}',
    stats.get('Number of reads', 'NA'),
    stats.get('Yield [Gb]', 'NA'),
    stats.get('Yield [Gb] (>25kb)', 'NA'),
    stats.get('N50', 'NA'),
    stats.get('Median length', 'NA'),
    stats.get('Mean length', 'NA'),
    stats.get('Median identity', 'NA'),
]
print('\\t'.join(row))
" >> ${prefix}_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mqc.tsv
    """
}
