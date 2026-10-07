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

def open_maybe_gzip(path):
    with open(path, 'rb') as probe:
        magic = probe.read(2)
    return gzip.open(path, 'rt') if magic == bytes([0x1f, 0x8b]) else open(path, 'rt')

def cov_at_cum_fraction(cum_fraction_by_cov, fraction):
    # Highest coverage at which at least `fraction` of bases are still covered.
    for cov, frac in sorted(cum_fraction_by_cov.items(), reverse=True):
        if frac >= fraction:
            return cov
    return None

with open_maybe_gzip('${global_dist}') as fh:
    cum_fraction_by_cov = {}
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) == 3 and parts[0] == 'total':
            cum_fraction_by_cov[int(parts[1])] = float(parts[2])

median = cov_at_cum_fraction(cum_fraction_by_cov, 0.5)
q1     = cov_at_cum_fraction(cum_fraction_by_cov, 0.75)
q3     = cov_at_cum_fraction(cum_fraction_by_cov, 0.25)

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
    "#     - coverage_iqr_cv:",
    "#         title: 'IQR/Median'",
    "#         description: 'Coverage uniformity: interquartile range of the coverage distribution divided by median coverage. Lower is more uniform.'",
    "#         min: 0",
    "#         scale: 'BuPu'",
    "#         format: '{:,.2f}'",
    "#         hidden: true",
    "#         namespace: 'Mosdepth'",
    "Sample\\tiqr_coverage\\tcoverage_iqr_cv",
    "${meta.id}\\t" + iqr_str + "\\t" + iqr_cv_str,
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
