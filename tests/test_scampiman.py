import sys
import pathlib
import logging
import pytest

TEST_DATA = pathlib.Path(__file__).parent / "test_data"

EXPECTED_STATS = {
    "total_reads":            10,
    "unmapped":               2,
    "removed_reads_primary":  3,
    "kept_primary":           5,
    "kept_secondary":         0,
    "kept_supplementary":     2,
}


@pytest.fixture(autouse=True)
def _clean_logger():
    """Clear pct_logger handlers after each test to prevent accumulation."""
    yield
    logging.getLogger("pct_logger").handlers.clear()


def test_single_end_bam(tmp_path, monkeypatch):
    """Integration test: single-end BAM → full scampiman pipeline outputs."""
    monkeypatch.chdir(tmp_path)  # log file lands in tmp_path, not the repo root

    out_dir  = tmp_path / "results"
    temp_dir = tmp_path / "tmp"
    sample   = "test"

    monkeypatch.setattr(sys, "argv", [
        "scampiman",
        "-r",  str(TEST_DATA / "test.bam"),
        "-g",  str(TEST_DATA / "sars-cov-2_NC_045512.2.fasta"),
        "-b",  str(TEST_DATA / "sars-cov-2_ARTICv5.3.2_subset.primer.bed"),
        "-f",  "bam",
        "-t",  "files",
        "-c",  "single-end",
        "--temp", str(temp_dir),
        "-o",  str(out_dir),
        "-s",  sample,
        "--keep", "all",
    ])

    from scampiman.scampiman import scampiman
    scampiman()

    # All four TSV outputs must exist and have more than one line
    tsv_files = [
        out_dir / f"{sample}.summarystats.tsv",
        out_dir / f"{sample}.ampliconstats.tsv",
        out_dir / f"{sample}.amplicontable.tsv",
        out_dir / f"{sample}.samcov.tsv",
    ]
    for tsv in tsv_files:
        assert tsv.exists(), f"Missing output: {tsv.name}"
        lines = [l for l in tsv.read_text().splitlines() if l.strip()]
        assert len(lines) > 1, f"{tsv.name} has only {len(lines)} line(s)"

    # summarystats values must match expected flagstats exactly
    lines = [l for l in (out_dir / f"{sample}.summarystats.tsv").read_text().splitlines() if l.strip()]
    header = lines[0].split("\t")
    values = [int(v) for v in lines[1].split("\t")]
    stats  = dict(zip(header, values))

    for key, expected in EXPECTED_STATS.items():
        assert stats[key] == expected, f"summarystats[{key}]: expected {expected}, got {stats.get(key)}"
