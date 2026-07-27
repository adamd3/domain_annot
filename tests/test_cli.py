"""Unit tests for CLI interface."""

from click.testing import CliRunner
from domain_annot.cli import main


def test_cli_version():
    runner = CliRunner()
    result = runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert "domain-annot" in result.output
    assert "0.1.0" in result.output


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "fetch" in result.output
    assert "process" in result.output


def test_process_help():
    runner = CliRunner()
    result = runner.invoke(main, ["process", "--help"])
    assert result.exit_code == 0
    assert "-o, --output" in result.output
    assert "Output path prefix" in result.output


def test_process_cmd_output_path():
    runner = CliRunner()
    with runner.isolated_filesystem():
        # Write dummy interpro tsv & entry list
        with open("interpro.tsv", "w") as f:
            f.write("prot1\thash\t100\tPfam\tPF00001\tdesc\t1\t50\t1e-5\tT\t01-01-2025\tIPR000001\tEntryDesc\t1\t50\n")
        with open("entry_list.txt", "w") as f:
            f.write("IPR000001\tDomain\tEntryDesc\n")

        result = runner.invoke(main, [
            "process",
            "-i", "interpro.tsv",
            "-e", "entry_list.txt",
            "-o", "custom_dir/my_output"
        ])
        assert result.exit_code == 0
        assert "Saved TSV domain annotations: custom_dir/my_output.tsv" in result.output


