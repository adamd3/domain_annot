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
