# scEnhancers

Generic Python workflow for single-cell enhancer discovery:

1. Read input metadata/FASTQ locations.
2. Align to reference genome with `cellranger count`.
3. Read enhancer-region output and retain high-confidence regions.

## Minimal usage

```python
from pathlib import Path
from scenhancers import ScEnhancerWorkflow, WorkflowConfig

config = WorkflowConfig(
    fastq_dir=Path("/data/fastq"),
    sample_id="sample-1",
    transcriptome=Path("/refs/grch38"),
    run_id="run-1",
    output_dir=Path("/output"),
    min_enhancer_score=0.5,
)

workflow = ScEnhancerWorkflow(config)

# runner can wrap subprocess.run for real execution
regions = workflow.run(Path("/path/to/enhancers.tsv"), runner=lambda cmd: print(cmd))
```

Enhancer input files must include these columns: `chromosome`, `start`, `end`, `score`.
TSV and CSV are both supported.
