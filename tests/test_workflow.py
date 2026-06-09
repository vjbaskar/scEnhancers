from pathlib import Path
import tempfile
import unittest

from scenhancers import ScEnhancerWorkflow, WorkflowConfig


class WorkflowTests(unittest.TestCase):
    def _config(self) -> WorkflowConfig:
        return WorkflowConfig(
            fastq_dir=Path("/data/fastq"),
            sample_id="sample-1",
            transcriptome=Path("/refs/grch38"),
            run_id="run-1",
            output_dir=Path("/output"),
            min_enhancer_score=0.5,
        )

    def test_cellranger_command_includes_required_args(self) -> None:
        workflow = ScEnhancerWorkflow(self._config())

        command = workflow.build_cellranger_command()

        self.assertEqual(command[0:2], ["cellranger", "count"])
        self.assertIn("--id=run-1", command)
        self.assertIn("--sample=sample-1", command)
        self.assertIn("--fastqs=/data/fastq", command)

    def test_read_enhancer_regions_filters_by_score(self) -> None:
        workflow = ScEnhancerWorkflow(self._config())

        with tempfile.TemporaryDirectory() as tmp:
            enhancer_file = Path(tmp) / "enhancers.tsv"
            enhancer_file.write_text(
                "chromosome\tstart\tend\tscore\n"
                "chr1\t100\t200\t0.9\n"
                "chr1\t300\t400\t0.1\n",
                encoding="utf-8",
            )

            regions = workflow.read_enhancer_regions(enhancer_file)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].chromosome, "chr1")
        self.assertEqual(regions[0].start, 100)

    def test_run_executes_alignment_then_reads_regions(self) -> None:
        workflow = ScEnhancerWorkflow(self._config())
        captured_command = []

        def mock_runner(command):
            captured_command.extend(command)

        with tempfile.TemporaryDirectory() as tmp:
            enhancer_file = Path(tmp) / "enhancers.csv"
            enhancer_file.write_text(
                "chromosome,start,end,score\n"
                "chr2,500,800,0.75\n",
                encoding="utf-8",
            )

            regions = workflow.run(enhancer_file, mock_runner)

        self.assertTrue(captured_command)
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].chromosome, "chr2")

    def test_read_enhancer_regions_raises_for_unparsable_score(self) -> None:
        workflow = ScEnhancerWorkflow(self._config())

        with tempfile.TemporaryDirectory() as tmp:
            enhancer_file = Path(tmp) / "enhancers.tsv"
            enhancer_file.write_text(
                "chromosome\tstart\tend\tscore\n"
                "chr1\t100\t200\tnot-a-number\n",
                encoding="utf-8",
            )

            with self.assertRaises(ValueError):
                workflow.read_enhancer_regions(enhancer_file)

    def test_read_enhancer_regions_raises_for_unknown_format(self) -> None:
        workflow = ScEnhancerWorkflow(self._config())

        with tempfile.TemporaryDirectory() as tmp:
            enhancer_file = Path(tmp) / "enhancers.txt"
            enhancer_file.write_text(
                "just some free-form text without delimiter hints",
                encoding="utf-8",
            )

            with self.assertRaises(ValueError):
                workflow.read_enhancer_regions(enhancer_file)


if __name__ == "__main__":
    unittest.main()
