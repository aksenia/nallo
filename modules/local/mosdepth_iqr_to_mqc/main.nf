process MOSDEPTH_IQR_TO_MQC {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(global_dist)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
python3 << 'PYEOF'
import gzip

with gzip.open('${global_dist}', 'rt') as fh:
    cum_fraction_by_cov = {}
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) == 3 and parts[0] == 'total':
            cum_fraction_by_cov[int(parts[1])] = float(parts[2])

median, q1, q3 = None, None, None
for cov, frac in sorted(cum_fraction_by_cov.items(), reverse=True):
    if median is None and frac >= 0.5:
        median = cov
    if q3 is None and frac >= 0.25:
        q3 = cov
    if q1 is None and frac >= 0.75:
        q1 = cov
        break

iqr    = (q3 - q1)      if q1 is not None and q3 is not None            else None
iqr_cv = (iqr / median) if iqr is not None and median not in (None, 0)  else None

iqr_str    = f"{iqr:.4f}"    if iqr    is not None else ''
iqr_cv_str = f"{iqr_cv:.4f}" if iqr_cv is not None else ''

lines = [
    "# id: 'mosdepth-iqr-coverage'",
    "# plot_type: 'generalstats'",
    "# pconfig:",
    "#     - iqr_coverage:",
    "#         title: 'IQR'",
    "#         description: 'Interquartile range of coverage depth (Q3 - Q1)'",
    "#         min: 0",
    "#         suffix: 'X'",
    "#         scale: 'BuPu'",
    "#         namespace: 'Mosdepth'",
    "#     - coefficient_of_iqr_variance:",
    "#         title: 'IQR CV'",
    "#         description: 'Coefficient of IQR variance (IQR / median). Lower = more uniform coverage.'",
    "#         min: 0",
    "#         scale: 'RdYlGn-rev'",
    "#         namespace: 'Mosdepth'",
    "Sample\\tiqr_coverage\\tcoefficient_of_iqr_variance",
    "${prefix}\\t" + iqr_str + "\\t" + iqr_cv_str,
]
with open('${prefix}_mqc.tsv', 'w') as fh:
    fh.write('\\n'.join(lines) + '\\n')
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mqc.tsv
    """
}
