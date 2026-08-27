process MOSDEPTH_REGIONS_TO_MQC {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(regions_bed), path(thresholds_bed)

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
from collections import Counter

def parse_regions(path):
    by_region = {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            fields = line.rstrip('\\n').split('\\t')
            if not fields or fields[0].startswith('#'):
                continue
            if len(fields) == 5:
                chrom, start, end, name, mean = fields
                by_region[(chrom, int(start), int(end))] = (name, float(mean))
            elif len(fields) == 4:
                chrom, start, end, mean = fields
                by_region[(chrom, int(start), int(end))] = (None, float(mean))
    return by_region

def parse_thresholds(path):
    with gzip.open(path, 'rt') as fh:
        lines = [l.rstrip('\\n') for l in fh]
    if not lines:
        return [], {}
    header = lines[0].lstrip('#').split('\\t')
    threshold_vals = [int(col.rstrip('X')) for col in header[4:]]
    by_region = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split('\\t')
        chrom, start, end, name = fields[:4]
        counts = [int(x) for x in fields[4:]]
        by_region[(chrom, int(start), int(end))] = (name, counts)
    return threshold_vals, by_region

sample_name = '${prefix}'
regions_data = parse_regions('${regions_bed}')
threshold_vals, thresh_data = parse_thresholds('${thresholds_bed}')

region_keys = list(thresh_data.keys())
for key in regions_data:
    if key not in thresh_data:
        region_keys.append(key)

def get_region_name(key):
    t_name = thresh_data[key][0] if key in thresh_data else None
    r_name = regions_data[key][0] if key in regions_data else None
    name = t_name or r_name
    return name if name is not None else f"{key[0]}:{key[1]}-{key[2]}"

name_counts = Counter(get_region_name(k) for k in region_keys)
duplicate_names = {name for name, count in name_counts.items() if count > 1}

out_lines = [
    "# id: 'mosdepth-per-region-coverage'",
    "# section_name: 'Per-region coverage'",
    "# description: 'Mean coverage and the percentage of bases at each requested threshold, per target region. Generated when mosdepth is run with --by BED_FILE and --thresholds.'",
    "# plot_type: 'table'",
    "# pconfig:",
    "#     id: 'mosdepth-per-region-table'",
    "#     title: 'Mosdepth: Per-region coverage'",
    "#     col1_header: 'Sample | Region'",
    "# headers:",
    "#     coordinates:",
    "#         title: 'Coordinates'",
    "#         description: 'Chromosome:start-end'",
    "#         scale: False",
    "#     mean_coverage:",
    "#         title: 'Mean Cov.'",
    "#         description: 'Mean coverage across the region'",
    "#         min: 0",
    "#         suffix: 'X'",
    "#         scale: 'BuPu'",
]
for t in threshold_vals:
    out_lines += [
        f"#     pct_at_{t}x:",
        f"#         title: '\\u2265 {t}X'",
        f"#         description: 'Percentage of bases in the region covered at least {t}X'",
        "#         min: 0",
        "#         max: 100",
        "#         suffix: '%'",
        "#         scale: 'RdYlGn'",
    ]

col_keys = ['coordinates', 'mean_coverage'] + [f'pct_at_{t}x' for t in threshold_vals]
out_lines.append('Sample | Region\\t' + '\\t'.join(col_keys))

for key in region_keys:
    chrom, start, end = key
    name = get_region_name(key)
    display_name = f"{name} ({chrom}:{start}-{end})" if name in duplicate_names else name
    row_key = f"{sample_name} | {display_name}"
    coords = f"{chrom}:{start}-{end}"
    mean = f"{regions_data[key][1]:.2f}" if key in regions_data else ''
    if key in thresh_data:
        length = end - start
        _, counts = thresh_data[key]
        pct_cols = [str(round(100.0 * c / length, 2) if length > 0 else 0.0) for c in counts]
    else:
        pct_cols = [''] * len(threshold_vals)
    out_lines.append('\\t'.join([row_key, coords, mean] + pct_cols))

with open('${prefix}_mqc.tsv', 'w') as fh:
    fh.write('\\n'.join(out_lines) + '\\n')
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mqc.tsv
    """
}
