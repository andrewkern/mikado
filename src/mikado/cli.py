"""Command-line interface for Mikado."""

from __future__ import annotations

from pathlib import Path

import click
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.style import Style

from mikado import __version__
from mikado.analysis.asymptotic import asymptotic_mk_test
from mikado.analysis.mk_test import mk_test
from mikado.analysis.polarized import polarized_mk_test
from mikado.io.output import OutputFormat, format_batch_results, format_result


class RainbowBarColumn(BarColumn):
    """A progress bar that cycles through rainbow colors."""

    RAINBOW_COLORS = [
        "#FF0000",  # Red
        "#FF7F00",  # Orange
        "#FFFF00",  # Yellow
        "#00FF00",  # Green
        "#0000FF",  # Blue
        "#4B0082",  # Indigo
        "#9400D3",  # Violet
    ]

    def __init__(self) -> None:
        super().__init__(bar_width=40)
        self._color_index = 0

    def render(self, task):  # type: ignore[no-untyped-def]
        """Render the bar with rainbow colors."""
        # Cycle through colors based on completed percentage
        if task.total:
            progress = task.completed / task.total
            color_idx = int(progress * len(self.RAINBOW_COLORS) * 3) % len(
                self.RAINBOW_COLORS
            )
        else:
            color_idx = self._color_index
            self._color_index = (self._color_index + 1) % len(self.RAINBOW_COLORS)

        self.complete_style = Style(color=self.RAINBOW_COLORS[color_idx])
        self.finished_style = Style(color="#9400D3")  # Violet when done
        return super().render(task)


def create_rainbow_progress() -> Progress:
    """Create a rainbow-colored progress bar."""
    return Progress(
        SpinnerColumn(style="bold magenta"),
        TextColumn("[bold blue]{task.description}"),
        RainbowBarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
    )


@click.group()
@click.version_option(version=__version__)
def main() -> None:
    """Mikado: McDonald-Kreitman test toolkit.

    A modern Python implementation for detecting selection using
    the McDonald-Kreitman test and related methods.
    """
    pass


@main.command()
@click.argument("ingroup", type=click.Path(exists=True, path_type=Path))
@click.argument("outgroup", type=click.Path(exists=True, path_type=Path), required=False)
@click.option(
    "-p",
    "--polarize",
    type=click.Path(exists=True, path_type=Path),
    help="Second outgroup file for polarized MK test",
)
@click.option(
    "--ingroup-match",
    help="Filter pattern for ingroup sequences (combined file mode)",
)
@click.option(
    "--outgroup-match",
    help="Filter pattern for outgroup sequences (combined file mode)",
)
@click.option(
    "--polarize-match",
    help="Filter pattern for second outgroup (combined file polarized mode)",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["pretty", "tsv", "json"]),
    default="pretty",
    help="Output format",
)
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
@click.option(
    "--pool-polymorphisms",
    is_flag=True,
    help="Pool polymorphisms from both populations (libsequence convention)",
)
def test(
    ingroup: Path,
    outgroup: Path | None,
    polarize: Path | None,
    ingroup_match: str | None,
    outgroup_match: str | None,
    polarize_match: str | None,
    output_format: str,
    reading_frame: int,
    pool_polymorphisms: bool,
) -> None:
    """Run the McDonald-Kreitman test.

    INGROUP: Path to FASTA file (ingroup sequences, or combined alignment)

    OUTGROUP: Path to FASTA file with outgroup sequences (optional if using combined mode)

    Examples:

        mikado test ingroup.fa outgroup.fa

        mikado test ingroup.fa outgroup.fa -p outgroup2.fa

        mikado test combined.fa --ingroup-match "gamb" --outgroup-match "002019"

        mikado test combined.fa --ingroup-match "gamb" --outgroup-match "002019" --polarize-match "amin"
    """
    from mikado.core.sequences import SequenceSet

    fmt = OutputFormat(output_format)

    # Determine mode: separate files or combined file
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined file mode
        if not ingroup_match or not outgroup_match:
            click.echo(
                "Error: Combined mode requires both --ingroup-match and --outgroup-match",
                err=True,
            )
            raise SystemExit(1)

        if outgroup is not None:
            click.echo(
                "Error: Do not provide OUTGROUP file when using --ingroup-match/--outgroup-match",
                err=True,
            )
            raise SystemExit(1)

        # Load combined file and filter
        all_seqs = SequenceSet.from_fasta(ingroup, reading_frame=reading_frame)
        ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
        outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

        if len(ingroup_seqs) == 0:
            click.echo(
                f"Error: No sequences match ingroup pattern '{ingroup_match}'", err=True
            )
            raise SystemExit(1)
        if len(outgroup_seqs) == 0:
            click.echo(
                f"Error: No sequences match outgroup pattern '{outgroup_match}'", err=True
            )
            raise SystemExit(1)

        click.echo(
            f"Combined mode: {len(ingroup_seqs)} ingroup, {len(outgroup_seqs)} outgroup sequences",
            err=True,
        )

        if polarize_match:
            outgroup2_seqs = all_seqs.filter_by_name(polarize_match)
            if len(outgroup2_seqs) == 0:
                click.echo(
                    f"Error: No sequences match polarize pattern '{polarize_match}'", err=True
                )
                raise SystemExit(1)
            click.echo(f"Polarizing with {len(outgroup2_seqs)} outgroup2 sequences", err=True)
            result = polarized_mk_test(
                ingroup=ingroup_seqs,
                outgroup1=outgroup_seqs,
                outgroup2=outgroup2_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
            )
        else:
            result = mk_test(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
            )
    else:
        # Separate files mode (original behavior)
        if outgroup is None:
            click.echo(
                "Error: OUTGROUP file required (or use --ingroup-match/--outgroup-match for combined mode)",
                err=True,
            )
            raise SystemExit(1)

        if polarize:
            result = polarized_mk_test(
                ingroup=ingroup,
                outgroup1=outgroup,
                outgroup2=polarize,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
            )
        else:
            result = mk_test(
                ingroup=ingroup,
                outgroup=outgroup,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
            )

    click.echo(format_result(result, fmt))


@main.command()
@click.argument("ingroup", type=click.Path(exists=True, path_type=Path))
@click.argument("outgroup", type=click.Path(exists=True, path_type=Path), required=False)
@click.option(
    "--ingroup-match",
    help="Filter pattern for ingroup sequences (combined file mode)",
)
@click.option(
    "--outgroup-match",
    help="Filter pattern for outgroup sequences (combined file mode)",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["pretty", "tsv", "json"]),
    default="pretty",
    help="Output format",
)
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
@click.option(
    "-b",
    "--bins",
    type=int,
    default=10,
    help="Number of frequency bins",
)
@click.option(
    "--bootstrap",
    type=int,
    default=100,
    help="Number of bootstrap replicates for CI",
)
@click.option(
    "--pool-polymorphisms",
    is_flag=True,
    help="Pool polymorphisms from both populations (libsequence convention)",
)
def asymptotic(
    ingroup: Path,
    outgroup: Path | None,
    ingroup_match: str | None,
    outgroup_match: str | None,
    output_format: str,
    reading_frame: int,
    bins: int,
    bootstrap: int,
    pool_polymorphisms: bool,
) -> None:
    """Run the asymptotic McDonald-Kreitman test.

    This method accounts for slightly deleterious mutations by examining
    how alpha changes across the frequency spectrum.

    Based on Messer & Petrov (2013) PNAS.

    INGROUP: Path to FASTA file (ingroup sequences, or combined alignment)

    OUTGROUP: Path to FASTA file with outgroup sequences (optional if using combined mode)

    Examples:

        mikado asymptotic ingroup.fa outgroup.fa

        mikado asymptotic ingroup.fa outgroup.fa --bins 20

        mikado asymptotic combined.fa --ingroup-match "gamb" --outgroup-match "002019"
    """
    from mikado.core.sequences import SequenceSet

    fmt = OutputFormat(output_format)

    # Determine mode: separate files or combined file
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined file mode
        if not ingroup_match or not outgroup_match:
            click.echo(
                "Error: Combined mode requires both --ingroup-match and --outgroup-match",
                err=True,
            )
            raise SystemExit(1)

        if outgroup is not None:
            click.echo(
                "Error: Do not provide OUTGROUP file when using --ingroup-match/--outgroup-match",
                err=True,
            )
            raise SystemExit(1)

        # Load combined file and filter
        all_seqs = SequenceSet.from_fasta(ingroup, reading_frame=reading_frame)
        ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
        outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

        if len(ingroup_seqs) == 0:
            click.echo(
                f"Error: No sequences match ingroup pattern '{ingroup_match}'", err=True
            )
            raise SystemExit(1)
        if len(outgroup_seqs) == 0:
            click.echo(
                f"Error: No sequences match outgroup pattern '{outgroup_match}'", err=True
            )
            raise SystemExit(1)

        click.echo(
            f"Combined mode: {len(ingroup_seqs)} ingroup, {len(outgroup_seqs)} outgroup sequences",
            err=True,
        )

        result = asymptotic_mk_test(
            ingroup=ingroup_seqs,
            outgroup=outgroup_seqs,
            reading_frame=reading_frame,
            num_bins=bins,
            bootstrap_replicates=bootstrap,
            pool_polymorphisms=pool_polymorphisms,
        )
    else:
        # Separate files mode (original behavior)
        if outgroup is None:
            click.echo(
                "Error: OUTGROUP file required (or use --ingroup-match/--outgroup-match for combined mode)",
                err=True,
            )
            raise SystemExit(1)

        result = asymptotic_mk_test(
            ingroup=ingroup,
            outgroup=outgroup,
            reading_frame=reading_frame,
            num_bins=bins,
            bootstrap_replicates=bootstrap,
            pool_polymorphisms=pool_polymorphisms,
        )

    click.echo(format_result(result, fmt))


@main.command()
@click.argument("input_dir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--file-pattern",
    default="*.fa",
    help="Glob pattern to match alignment files (combined mode)",
)
@click.option(
    "--outgroup-pattern",
    default="*_outgroup.fa",
    help="Glob pattern to match outgroup files (separate files mode)",
)
@click.option(
    "--ingroup-pattern",
    default="*_ingroup.fa",
    help="Glob pattern to match ingroup files (separate files mode)",
)
@click.option(
    "--ingroup-match",
    help="Filter pattern for ingroup sequences (combined file mode)",
)
@click.option(
    "--outgroup-match",
    help="Filter pattern for outgroup sequences (combined file mode)",
)
@click.option(
    "--polarize-match",
    help="Filter pattern for second outgroup (combined file polarized mode)",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["pretty", "tsv", "json"]),
    default="tsv",
    help="Output format",
)
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
@click.option(
    "-a",
    "--asymptotic",
    is_flag=True,
    help="Use asymptotic MK test instead of standard",
)
@click.option(
    "-p",
    "--polarize-pattern",
    default=None,
    help="Glob pattern for second outgroup files (separate files polarized mode)",
)
@click.option(
    "-b",
    "--bins",
    type=int,
    default=10,
    help="Number of frequency bins (asymptotic only)",
)
@click.option(
    "--bootstrap",
    type=int,
    default=100,
    help="Bootstrap replicates for CI (asymptotic only)",
)
@click.option(
    "--pool-polymorphisms",
    is_flag=True,
    help="Pool polymorphisms from both populations (libsequence convention)",
)
def batch(
    input_dir: Path,
    file_pattern: str,
    outgroup_pattern: str,
    ingroup_pattern: str,
    ingroup_match: str | None,
    outgroup_match: str | None,
    polarize_match: str | None,
    output_format: str,
    reading_frame: int,
    asymptotic: bool,
    polarize_pattern: str | None,
    bins: int,
    bootstrap: int,
    pool_polymorphisms: bool,
) -> None:
    """Run MK test on multiple gene files.

    Supports two modes:

    1. Separate files mode (default): Pairs of ingroup/outgroup FASTA files

    2. Combined files mode: Single alignment files with sequences from multiple species,
       filtered by --ingroup-match and --outgroup-match patterns

    INPUT_DIR: Directory containing FASTA files

    Examples (separate files mode):

        mikado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

        mikado batch genes/ --asymptotic --bins 10

    Examples (combined files mode):

        mikado batch alignments/ --file-pattern "*.fasta" --ingroup-match "gamb" --outgroup-match "002019"

        mikado batch alignments/ --file-pattern "*.fasta" --ingroup-match "gamb" --outgroup-match "002019" --polarize-match "amin"
    """
    from mikado.core.sequences import SequenceSet

    fmt = OutputFormat(output_format)

    # Determine mode: combined files or separate files
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined files mode
        if not ingroup_match or not outgroup_match:
            click.echo(
                "Error: Combined mode requires both --ingroup-match and --outgroup-match",
                err=True,
            )
            return

        # Check for mutually exclusive options
        if asymptotic and polarize_match:
            click.echo(
                "Error: --asymptotic and --polarize-match are mutually exclusive", err=True
            )
            return

        # Find all alignment files
        alignment_files = sorted(input_dir.glob(file_pattern))

        if not alignment_files:
            click.echo(
                f"No files matching '{file_pattern}' found in {input_dir}", err=True
            )
            return

        results = []
        warnings = []

        with create_rainbow_progress() as progress:
            task = progress.add_task(
                "Processing alignments", total=len(alignment_files)
            )

            for alignment_file in alignment_files:
                try:
                    all_seqs = SequenceSet.from_fasta(
                        alignment_file, reading_frame=reading_frame
                    )
                    ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
                    outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

                    if len(ingroup_seqs) == 0:
                        warnings.append(
                            f"Warning: No ingroup sequences in {alignment_file.name}"
                        )
                        progress.advance(task)
                        continue
                    if len(outgroup_seqs) == 0:
                        warnings.append(
                            f"Warning: No outgroup sequences in {alignment_file.name}"
                        )
                        progress.advance(task)
                        continue

                    if asymptotic:
                        result = asymptotic_mk_test(
                            ingroup=ingroup_seqs,
                            outgroup=outgroup_seqs,
                            reading_frame=reading_frame,
                            num_bins=bins,
                            bootstrap_replicates=bootstrap,
                            pool_polymorphisms=pool_polymorphisms,
                        )
                    elif polarize_match:
                        outgroup2_seqs = all_seqs.filter_by_name(polarize_match)
                        if len(outgroup2_seqs) == 0:
                            warnings.append(
                                f"Warning: No outgroup2 sequences in {alignment_file.name}"
                            )
                            progress.advance(task)
                            continue
                        result = polarized_mk_test(
                            ingroup=ingroup_seqs,
                            outgroup1=outgroup_seqs,
                            outgroup2=outgroup2_seqs,
                            reading_frame=reading_frame,
                            pool_polymorphisms=pool_polymorphisms,
                        )
                    else:
                        result = mk_test(
                            ingroup=ingroup_seqs,
                            outgroup=outgroup_seqs,
                            reading_frame=reading_frame,
                            pool_polymorphisms=pool_polymorphisms,
                        )
                    results.append((alignment_file.stem, result))
                except Exception as e:
                    warnings.append(f"Error processing {alignment_file.name}: {e}")

                progress.advance(task)

        # Print warnings after progress bar completes
        for warning in warnings:
            click.echo(warning, err=True)

    else:
        # Separate files mode (original behavior)
        # Check for mutually exclusive options
        if asymptotic and polarize_pattern:
            click.echo(
                "Error: --asymptotic and --polarize-pattern are mutually exclusive", err=True
            )
            return

        # Find all ingroup files
        ingroup_files = sorted(input_dir.glob(ingroup_pattern))

        if not ingroup_files:
            click.echo(
                f"No files matching '{ingroup_pattern}' found in {input_dir}", err=True
            )
            return

        results = []
        warnings = []

        with create_rainbow_progress() as progress:
            task = progress.add_task("Processing files", total=len(ingroup_files))

            for ingroup_file in ingroup_files:
                # Construct outgroup filename
                base_name = ingroup_file.stem
                if "_ingroup" in base_name:
                    outgroup_name = base_name.replace("_ingroup", "_outgroup") + ".fa"
                elif "_in" in base_name:
                    outgroup_name = base_name.replace("_in", "_out") + ".fa"
                else:
                    outgroup_name = base_name + "_outgroup.fa"

                outgroup_file = input_dir / outgroup_name

                if not outgroup_file.exists():
                    potential_matches = list(input_dir.glob(outgroup_pattern))
                    prefix = base_name.split("_")[0]
                    matches = [
                        f for f in potential_matches if f.stem.startswith(prefix)
                    ]
                    if matches:
                        outgroup_file = matches[0]
                    else:
                        warnings.append(
                            f"Warning: No outgroup found for {ingroup_file.name}"
                        )
                        progress.advance(task)
                        continue

                outgroup2_file: Path | None = None
                if polarize_pattern:
                    potential_matches = list(input_dir.glob(polarize_pattern))
                    prefix = base_name.split("_")[0]
                    matches = [
                        f for f in potential_matches if f.stem.startswith(prefix)
                    ]
                    if matches:
                        outgroup2_file = matches[0]
                    else:
                        warnings.append(
                            f"Warning: No second outgroup found for {ingroup_file.name}"
                        )
                        progress.advance(task)
                        continue

                try:
                    if asymptotic:
                        result = asymptotic_mk_test(
                            ingroup=ingroup_file,
                            outgroup=outgroup_file,
                            reading_frame=reading_frame,
                            num_bins=bins,
                            bootstrap_replicates=bootstrap,
                            pool_polymorphisms=pool_polymorphisms,
                        )
                    elif polarize_pattern and outgroup2_file:
                        result = polarized_mk_test(
                            ingroup=ingroup_file,
                            outgroup1=outgroup_file,
                            outgroup2=outgroup2_file,
                            reading_frame=reading_frame,
                            pool_polymorphisms=pool_polymorphisms,
                        )
                    else:
                        result = mk_test(
                            ingroup=ingroup_file,
                            outgroup=outgroup_file,
                            reading_frame=reading_frame,
                            pool_polymorphisms=pool_polymorphisms,
                        )
                    results.append((ingroup_file.stem, result))
                except Exception as e:
                    warnings.append(f"Error processing {ingroup_file.name}: {e}")

                progress.advance(task)

        # Print warnings after progress bar completes
        for warning in warnings:
            click.echo(warning, err=True)

    if results:
        click.echo(format_batch_results(results, fmt))
    else:
        click.echo("No results to display", err=True)


@main.command()
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
def info(fasta: Path, reading_frame: int) -> None:
    """Display information about a FASTA file.

    Shows sequence count, alignment length, and codon statistics.

    FASTA: Path to FASTA file
    """
    from mikado.core.sequences import SequenceSet

    seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)

    click.echo(f"File: {fasta.name}")
    click.echo(f"Sequences: {len(seqs)}")

    if seqs.sequences:
        click.echo(f"Alignment length: {seqs.alignment_length} bp")
        click.echo(f"Codons: {seqs.num_codons}")
        click.echo(f"Reading frame: {reading_frame}")

        # Count polymorphic sites
        poly_sites = seqs.polymorphic_codons()
        click.echo(f"Polymorphic codons: {len(poly_sites)}")

        # List sequence names
        click.echo("\nSequences:")
        for seq in seqs.sequences:
            click.echo(f"  {seq.name} ({len(seq)} bp)")


if __name__ == "__main__":
    main()
