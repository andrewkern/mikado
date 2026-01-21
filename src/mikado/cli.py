"""Command-line interface for Mikado."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated, Optional

import typer
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
from mikado.analysis.asymptotic import (
    asymptotic_mk_test,
    asymptotic_mk_test_aggregated,
)
from mikado.analysis.mk_test import mk_test
from mikado.analysis.polarized import polarized_mk_test
from mikado.batch_workers import BatchTask, WorkerResult, process_gene
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


def get_worker_count(requested: int, num_tasks: int) -> int:
    """Determine the optimal number of workers.

    Args:
        requested: Number of workers requested (0=auto, 1=sequential)
        num_tasks: Number of tasks to process

    Returns:
        Number of workers to use
    """
    # Sequential mode if explicitly requested or too few tasks
    if requested == 1 or num_tasks < 10:
        return 1

    cpu_count = os.cpu_count() or 4

    if requested > 0:
        # User requested specific count, cap at CPU count
        return min(requested, cpu_count)

    # Auto mode: use cpu_count - 1 (leave one for main process)
    # but cap at number of tasks
    return max(1, min(cpu_count - 1, num_tasks))


def run_parallel_batch(
    tasks: list[BatchTask],
    num_workers: int,
    description: str,
) -> tuple[list[WorkerResult], list[str]]:
    """Run batch processing with ProcessPoolExecutor.

    Args:
        tasks: List of BatchTask objects to process
        num_workers: Number of worker processes (1 = sequential)
        description: Description for progress bar

    Returns:
        Tuple of (results, warnings) where results are WorkerResult objects
        and warnings are string messages
    """
    results: list[WorkerResult] = []
    warnings: list[str] = []

    if num_workers == 1:
        # Sequential mode - process in main process
        with create_rainbow_progress() as progress:
            task_id = progress.add_task(description, total=len(tasks))
            for task in tasks:
                worker_result = process_gene(task)
                if worker_result.error:
                    warnings.append(worker_result.error)
                elif worker_result.warning:
                    warnings.append(f"Warning: {worker_result.warning}")
                elif worker_result.result is not None:
                    results.append(worker_result)
                progress.advance(task_id)
    else:
        # Parallel mode
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_gene, task): task for task in tasks}

            with create_rainbow_progress() as progress:
                task_id = progress.add_task(description, total=len(tasks))

                for future in as_completed(futures):
                    try:
                        worker_result = future.result()
                        if worker_result.error:
                            warnings.append(worker_result.error)
                        elif worker_result.warning:
                            warnings.append(f"Warning: {worker_result.warning}")
                        elif worker_result.result is not None:
                            results.append(worker_result)
                    except Exception as e:
                        task = futures[future]
                        warnings.append(f"Error processing {task.file_path.name}: {e}")
                    progress.advance(task_id)

    return results, warnings


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        typer.echo(f"mikado {__version__}")
        raise typer.Exit()


app = typer.Typer(
    name="mikado",
    help="Mikado: McDonald-Kreitman test toolkit.\n\n"
    "A modern Python implementation for detecting selection using "
    "the McDonald-Kreitman test and related methods.",
    no_args_is_help=True,
)


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option("--version", "-v", callback=version_callback, is_eager=True),
    ] = False,
) -> None:
    """Mikado: McDonald-Kreitman test toolkit."""
    pass


@app.command()
def test(
    ingroup: Annotated[
        Path,
        typer.Argument(help="Path to FASTA file (ingroup sequences, or combined alignment)"),
    ],
    outgroup: Annotated[
        Optional[Path],
        typer.Argument(help="Path to FASTA file with outgroup sequences"),
    ] = None,
    polarize: Annotated[
        Optional[Path],
        typer.Option("--polarize", "-p", help="Second outgroup file for polarized MK test"),
    ] = None,
    ingroup_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for ingroup sequences (combined file mode)"),
    ] = None,
    outgroup_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for outgroup sequences (combined file mode)"),
    ] = None,
    polarize_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for second outgroup (combined file polarized mode)"),
    ] = None,
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format"),
    ] = "pretty",
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1, 2, or 3)"),
    ] = 1,
    pool_polymorphisms: Annotated[
        bool,
        typer.Option(help="Pool polymorphisms from both populations (libsequence convention)"),
    ] = False,
    min_freq: Annotated[
        float,
        typer.Option(
            "--min-freq",
            help="Minimum derived allele frequency for polymorphisms (0.0-1.0)",
        ),
    ] = 0.0,
) -> None:
    """Run the McDonald-Kreitman test.

    Examples:

        mikado test ingroup.fa outgroup.fa

        mikado test ingroup.fa outgroup.fa -p outgroup2.fa

        mikado test combined.fa --ingroup-match "gamb" --outgroup-match "002019"

        mikado test combined.fa --ingroup-match "gamb" --outgroup-match "002019" --polarize-match "amin"
    """
    from mikado.core.sequences import SequenceSet

    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'. Must be pretty, tsv, or json.", err=True)
        raise typer.Exit(1)

    fmt = OutputFormat(output_format)

    # Determine mode: separate files or combined file
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined file mode
        if not ingroup_match or not outgroup_match:
            typer.echo(
                "Error: Combined mode requires both --ingroup-match and --outgroup-match",
                err=True,
            )
            raise typer.Exit(1)

        if outgroup is not None:
            typer.echo(
                "Error: Do not provide OUTGROUP file when using --ingroup-match/--outgroup-match",
                err=True,
            )
            raise typer.Exit(1)

        # Load combined file and filter
        all_seqs = SequenceSet.from_fasta(ingroup, reading_frame=reading_frame)
        ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
        outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

        if len(ingroup_seqs) == 0:
            typer.echo(
                f"Error: No sequences match ingroup pattern '{ingroup_match}'", err=True
            )
            raise typer.Exit(1)
        if len(outgroup_seqs) == 0:
            typer.echo(
                f"Error: No sequences match outgroup pattern '{outgroup_match}'", err=True
            )
            raise typer.Exit(1)

        typer.echo(
            f"Combined mode: {len(ingroup_seqs)} ingroup, {len(outgroup_seqs)} outgroup sequences",
            err=True,
        )

        if polarize_match:
            outgroup2_seqs = all_seqs.filter_by_name(polarize_match)
            if len(outgroup2_seqs) == 0:
                typer.echo(
                    f"Error: No sequences match polarize pattern '{polarize_match}'", err=True
                )
                raise typer.Exit(1)
            typer.echo(f"Polarizing with {len(outgroup2_seqs)} outgroup2 sequences", err=True)
            result = polarized_mk_test(
                ingroup=ingroup_seqs,
                outgroup1=outgroup_seqs,
                outgroup2=outgroup2_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
            )
        else:
            result = mk_test(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
            )
    else:
        # Separate files mode (original behavior)
        if outgroup is None:
            typer.echo(
                "Error: OUTGROUP file required (or use --ingroup-match/--outgroup-match for combined mode)",
                err=True,
            )
            raise typer.Exit(1)

        if polarize:
            result = polarized_mk_test(
                ingroup=ingroup,
                outgroup1=outgroup,
                outgroup2=polarize,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
            )
        else:
            result = mk_test(
                ingroup=ingroup,
                outgroup=outgroup,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
            )

    typer.echo(format_result(result, fmt))


@app.command()
def asymptotic(
    ingroup: Annotated[
        Path,
        typer.Argument(help="Path to FASTA file (ingroup sequences, or combined alignment)"),
    ],
    outgroup: Annotated[
        Optional[Path],
        typer.Argument(help="Path to FASTA file with outgroup sequences"),
    ] = None,
    ingroup_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for ingroup sequences (combined file mode)"),
    ] = None,
    outgroup_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for outgroup sequences (combined file mode)"),
    ] = None,
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format"),
    ] = "pretty",
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1, 2, or 3)"),
    ] = 1,
    bins: Annotated[
        int,
        typer.Option("--bins", "-b", help="Number of frequency bins"),
    ] = 10,
    bootstrap: Annotated[
        int,
        typer.Option(help="Number of bootstrap replicates for CI"),
    ] = 100,
    pool_polymorphisms: Annotated[
        bool,
        typer.Option(help="Pool polymorphisms from both populations (libsequence convention)"),
    ] = False,
) -> None:
    """Run the asymptotic McDonald-Kreitman test.

    This method accounts for slightly deleterious mutations by examining
    how alpha changes across the frequency spectrum.

    Based on Messer & Petrov (2013) PNAS.

    Examples:

        mikado asymptotic ingroup.fa outgroup.fa

        mikado asymptotic ingroup.fa outgroup.fa --bins 20

        mikado asymptotic combined.fa --ingroup-match "gamb" --outgroup-match "002019"
    """
    from mikado.core.sequences import SequenceSet

    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'. Must be pretty, tsv, or json.", err=True)
        raise typer.Exit(1)

    fmt = OutputFormat(output_format)

    # Determine mode: separate files or combined file
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined file mode
        if not ingroup_match or not outgroup_match:
            typer.echo(
                "Error: Combined mode requires both --ingroup-match and --outgroup-match",
                err=True,
            )
            raise typer.Exit(1)

        if outgroup is not None:
            typer.echo(
                "Error: Do not provide OUTGROUP file when using --ingroup-match/--outgroup-match",
                err=True,
            )
            raise typer.Exit(1)

        # Load combined file and filter
        all_seqs = SequenceSet.from_fasta(ingroup, reading_frame=reading_frame)
        ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
        outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

        if len(ingroup_seqs) == 0:
            typer.echo(
                f"Error: No sequences match ingroup pattern '{ingroup_match}'", err=True
            )
            raise typer.Exit(1)
        if len(outgroup_seqs) == 0:
            typer.echo(
                f"Error: No sequences match outgroup pattern '{outgroup_match}'", err=True
            )
            raise typer.Exit(1)

        typer.echo(
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
            typer.echo(
                "Error: OUTGROUP file required (or use --ingroup-match/--outgroup-match for combined mode)",
                err=True,
            )
            raise typer.Exit(1)

        result = asymptotic_mk_test(
            ingroup=ingroup,
            outgroup=outgroup,
            reading_frame=reading_frame,
            num_bins=bins,
            bootstrap_replicates=bootstrap,
            pool_polymorphisms=pool_polymorphisms,
        )

    typer.echo(format_result(result, fmt))


@app.command()
def batch(
    input_dir: Annotated[
        Path,
        typer.Argument(help="Directory containing FASTA files"),
    ],
    file_pattern: Annotated[
        str,
        typer.Option(help="Glob pattern to match alignment files (combined mode)"),
    ] = "*.fa",
    outgroup_pattern: Annotated[
        str,
        typer.Option(help="Glob pattern to match outgroup files (separate files mode)"),
    ] = "*_outgroup.fa",
    ingroup_pattern: Annotated[
        str,
        typer.Option(help="Glob pattern to match ingroup files (separate files mode)"),
    ] = "*_ingroup.fa",
    ingroup_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for ingroup sequences (combined file mode)"),
    ] = None,
    outgroup_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for outgroup sequences (combined file mode)"),
    ] = None,
    polarize_match: Annotated[
        Optional[str],
        typer.Option(help="Filter pattern for second outgroup (combined file polarized mode)"),
    ] = None,
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format"),
    ] = "tsv",
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1, 2, or 3)"),
    ] = 1,
    use_asymptotic: Annotated[
        bool,
        typer.Option("--asymptotic", "-a", help="Use asymptotic MK test instead of standard"),
    ] = False,
    polarize_pattern: Annotated[
        Optional[str],
        typer.Option("--polarize-pattern", "-p", help="Glob pattern for second outgroup files (separate files polarized mode)"),
    ] = None,
    bins: Annotated[
        int,
        typer.Option("--bins", "-b", help="Number of frequency bins (asymptotic only)"),
    ] = 10,
    bootstrap: Annotated[
        int,
        typer.Option(help="Bootstrap replicates for CI (asymptotic only)"),
    ] = 100,
    pool_polymorphisms: Annotated[
        bool,
        typer.Option(help="Pool polymorphisms from both populations (libsequence convention)"),
    ] = False,
    aggregate: Annotated[
        bool,
        typer.Option(
            "--aggregate/--per-gene",
            help="Aggregate data across genes (default) vs per-gene analysis (asymptotic only)",
        ),
    ] = True,
    freq_cutoffs: Annotated[
        str,
        typer.Option(
            "--freq-cutoffs",
            help="Frequency range for fitting as 'low,high' (asymptotic only)",
        ),
    ] = "0.1,0.9",
    min_freq: Annotated[
        float,
        typer.Option(
            "--min-freq",
            help="Minimum derived allele frequency for polymorphisms (standard MK only)",
        ),
    ] = 0.0,
    workers: Annotated[
        int,
        typer.Option(
            "--workers",
            "-w",
            help="Number of parallel workers (0=auto, 1=sequential)",
            min=0,
        ),
    ] = 0,
) -> None:
    """Run MK test on multiple gene files.

    Supports two modes:

    1. Separate files mode (default): Pairs of ingroup/outgroup FASTA files

    2. Combined files mode: Single alignment files with sequences from multiple species,
       filtered by --ingroup-match and --outgroup-match patterns

    Examples (separate files mode):

        mikado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

        mikado batch genes/ --asymptotic --bins 10

    Examples (combined files mode):

        mikado batch alignments/ --file-pattern "*.fasta" --ingroup-match "gamb" --outgroup-match "002019"

        mikado batch alignments/ --file-pattern "*.fasta" --ingroup-match "gamb" --outgroup-match "002019" --polarize-match "amin"
    """
    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'. Must be pretty, tsv, or json.", err=True)
        raise typer.Exit(1)

    fmt = OutputFormat(output_format)

    # Parse frequency cutoffs
    try:
        cutoff_parts = freq_cutoffs.split(",")
        frequency_cutoffs = (float(cutoff_parts[0]), float(cutoff_parts[1]))
    except (ValueError, IndexError):
        typer.echo(
            f"Error: Invalid frequency cutoffs '{freq_cutoffs}'. Use format 'low,high'",
            err=True,
        )
        raise typer.Exit(1)

    # Determine mode: combined files or separate files
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined files mode
        if not ingroup_match or not outgroup_match:
            typer.echo(
                "Error: Combined mode requires both --ingroup-match and --outgroup-match",
                err=True,
            )
            return

        # Check for mutually exclusive options
        if use_asymptotic and polarize_match:
            typer.echo(
                "Error: --asymptotic and --polarize-match are mutually exclusive", err=True
            )
            return

        # Find all alignment files
        alignment_files = sorted(input_dir.glob(file_pattern))

        if not alignment_files:
            typer.echo(
                f"No files matching '{file_pattern}' found in {input_dir}", err=True
            )
            return

        # Determine number of workers
        num_workers = get_worker_count(workers, len(alignment_files))

        # Build task list
        tasks = [
            BatchTask(
                file_path=f,
                ingroup_match=ingroup_match,
                outgroup_match=outgroup_match,
                polarize_match=polarize_match,
                reading_frame=reading_frame,
                use_asymptotic=use_asymptotic,
                bins=bins,
                bootstrap=bootstrap,
                pool_polymorphisms=pool_polymorphisms,
                min_freq=min_freq,
                extract_only=(use_asymptotic and aggregate),
            )
            for f in alignment_files
        ]

        # Aggregated asymptotic mode: extract data, then aggregate
        if use_asymptotic and aggregate:
            from mikado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            # Print warnings
            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                # Run aggregated analysis
                result = asymptotic_mk_test_aggregated(
                    gene_data=gene_data_list,
                    num_bins=bins,
                    ci_replicates=bootstrap * 100,  # More replicates for Monte Carlo
                    frequency_cutoffs=frequency_cutoffs,
                )
                typer.echo(format_result(result, fmt))
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Per-gene mode
        worker_results, warnings = run_parallel_batch(
            tasks, num_workers, "Processing alignments"
        )

        # Print warnings after progress bar completes
        for warning in warnings:
            typer.echo(warning, err=True)

        # Convert to results format
        results = [(r.gene_id, r.result) for r in worker_results]

    else:
        # Separate files mode (original behavior)
        # Check for mutually exclusive options
        if use_asymptotic and polarize_pattern:
            typer.echo(
                "Error: --asymptotic and --polarize-pattern are mutually exclusive", err=True
            )
            return

        # Find all ingroup files
        ingroup_files = sorted(input_dir.glob(ingroup_pattern))

        if not ingroup_files:
            typer.echo(
                f"No files matching '{ingroup_pattern}' found in {input_dir}", err=True
            )
            return

        # Helper function to find outgroup file
        def find_outgroup_file(ingroup_file: Path) -> Path | None:
            base_name = ingroup_file.stem
            if "_ingroup" in base_name:
                outgroup_name = base_name.replace("_ingroup", "_outgroup") + ".fa"
            elif "_in" in base_name:
                outgroup_name = base_name.replace("_in", "_out") + ".fa"
            else:
                outgroup_name = base_name + "_outgroup.fa"

            outgroup_file_path = input_dir / outgroup_name

            if not outgroup_file_path.exists():
                potential_matches = list(input_dir.glob(outgroup_pattern))
                prefix = base_name.split("_")[0]
                matches = [
                    f for f in potential_matches if f.stem.startswith(prefix)
                ]
                if matches:
                    return matches[0]
                return None
            return outgroup_file_path

        # Helper function to find second outgroup file
        def find_outgroup2_file(ingroup_file: Path) -> Path | None:
            if not polarize_pattern:
                return None
            base_name = ingroup_file.stem
            potential_matches = list(input_dir.glob(polarize_pattern))
            prefix = base_name.split("_")[0]
            matches = [
                f for f in potential_matches if f.stem.startswith(prefix)
            ]
            if matches:
                return matches[0]
            return None

        # Determine number of workers
        num_workers = get_worker_count(workers, len(ingroup_files))

        # Build task list
        tasks = []
        pre_warnings = []
        for ingroup_file in ingroup_files:
            outgroup_file = find_outgroup_file(ingroup_file)
            if outgroup_file is None:
                pre_warnings.append(
                    f"Warning: No outgroup found for {ingroup_file.name}"
                )
                continue

            outgroup2_file: Path | None = None
            if polarize_pattern:
                outgroup2_file = find_outgroup2_file(ingroup_file)
                if outgroup2_file is None:
                    pre_warnings.append(
                        f"Warning: No second outgroup found for {ingroup_file.name}"
                    )
                    continue

            tasks.append(
                BatchTask(
                    file_path=ingroup_file,
                    outgroup_file=outgroup_file,
                    outgroup2_file=outgroup2_file,
                    reading_frame=reading_frame,
                    use_asymptotic=use_asymptotic,
                    bins=bins,
                    bootstrap=bootstrap,
                    pool_polymorphisms=pool_polymorphisms,
                    min_freq=min_freq,
                    extract_only=(use_asymptotic and aggregate),
                )
            )

        # Print pre-processing warnings
        for warning in pre_warnings:
            typer.echo(warning, err=True)

        if not tasks:
            typer.echo("No valid file pairs found", err=True)
            return

        # Aggregated asymptotic mode: extract data, then aggregate
        if use_asymptotic and aggregate:
            from mikado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            # Print warnings
            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                result = asymptotic_mk_test_aggregated(
                    gene_data=gene_data_list,
                    num_bins=bins,
                    ci_replicates=bootstrap * 100,
                    frequency_cutoffs=frequency_cutoffs,
                )
                typer.echo(format_result(result, fmt))
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Per-gene mode
        worker_results, warnings = run_parallel_batch(
            tasks, num_workers, "Processing files"
        )

        # Print warnings after progress bar completes
        for warning in warnings:
            typer.echo(warning, err=True)

        # Convert to results format
        results = [(r.gene_id, r.result) for r in worker_results]

    if results:
        typer.echo(format_batch_results(results, fmt))
    else:
        typer.echo("No results to display", err=True)


@app.command()
def info(
    fasta: Annotated[
        Path,
        typer.Argument(help="Path to FASTA file"),
    ],
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1, 2, or 3)"),
    ] = 1,
) -> None:
    """Display information about a FASTA file.

    Shows sequence count, alignment length, and codon statistics.
    """
    from mikado.core.sequences import SequenceSet

    seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)

    typer.echo(f"File: {fasta.name}")
    typer.echo(f"Sequences: {len(seqs)}")

    if seqs.sequences:
        typer.echo(f"Alignment length: {seqs.alignment_length} bp")
        typer.echo(f"Codons: {seqs.num_codons}")
        typer.echo(f"Reading frame: {reading_frame}")

        # Count polymorphic sites
        poly_sites = seqs.polymorphic_codons()
        typer.echo(f"Polymorphic codons: {len(poly_sites)}")

        # List sequence names
        typer.echo("\nSequences:")
        for seq in seqs.sequences:
            typer.echo(f"  {seq.name} ({len(seq)} bp)")


if __name__ == "__main__":
    app()
