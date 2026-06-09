from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Sequence
import csv


@dataclass(frozen=True)
class WorkflowConfig:
    fastq_dir: Path
    sample_id: str
    transcriptome: Path
    run_id: str
    output_dir: Path
    min_enhancer_score: float = 0.0


@dataclass(frozen=True)
class EnhancerRegion:
    chromosome: str
    start: int
    end: int
    score: float


class ScEnhancerWorkflow:
    """Generic workflow: ingest data, run cellranger, and detect enhancer regions."""

    def __init__(self, config: WorkflowConfig):
        self.config = config

    def build_cellranger_command(self) -> List[str]:
        out_path = self.config.output_dir / self.config.run_id
        return [
            "cellranger",
            "count",
            f"--id={self.config.run_id}",
            f"--transcriptome={self.config.transcriptome}",
            f"--fastqs={self.config.fastq_dir}",
            f"--sample={self.config.sample_id}",
            f"--output-dir={out_path}",
        ]

    def run_alignment(self, runner: Callable[[Sequence[str]], None]) -> None:
        runner(self.build_cellranger_command())

    def read_enhancer_regions(self, enhancer_file: Path) -> List[EnhancerRegion]:
        if not enhancer_file.exists():
            raise FileNotFoundError(f"Enhancer file not found: {enhancer_file}")

        regions: List[EnhancerRegion] = []
        try:
            with enhancer_file.open("r", newline="") as handle:
                sample = handle.read(1024)
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
                except csv.Error as exc:
                    raise ValueError(
                        f"Unable to detect CSV/TSV format for enhancer file: {enhancer_file}"
                    ) from exc

                handle.seek(0)
                reader = csv.DictReader(handle, dialect=dialect)
                required = {"chromosome", "start", "end", "score"}
                missing = required.difference(reader.fieldnames or ())
                if missing:
                    raise ValueError(
                        f"Enhancer file missing required columns: {', '.join(sorted(missing))}"
                    )

                first_data_line_number = 2
                for idx, row in enumerate(reader, start=first_data_line_number):
                    try:
                        score = float(row["score"])
                        start = int(row["start"])
                        end = int(row["end"])
                    except (TypeError, ValueError) as exc:
                        raise ValueError(
                            f"Failed parsing enhancer row {idx} in {enhancer_file}: {row}"
                        ) from exc

                    if score >= self.config.min_enhancer_score:
                        regions.append(
                            EnhancerRegion(
                                chromosome=row["chromosome"],
                                start=start,
                                end=end,
                                score=score,
                            )
                        )
        except OSError as exc:
            raise OSError(f"Unable to read enhancer file: {enhancer_file}") from exc

        return regions

    def run(
        self,
        enhancer_file: Path,
        runner: Callable[[Sequence[str]], None],
    ) -> List[EnhancerRegion]:
        self.run_alignment(runner)
        return self.read_enhancer_regions(enhancer_file)
