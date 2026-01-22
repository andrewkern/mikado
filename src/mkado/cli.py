"""Command-line interface for mkado."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.style import Style

from mkado import __version__
from mkado.analysis.asymptotic import (
    asymptotic_mk_test,
    asymptotic_mk_test_aggregated,
)
from mkado.analysis.mk_test import mk_test
from mkado.analysis.polarized import polarized_mk_test
from mkado.batch_workers import BatchTask, WorkerResult, process_gene
from mkado.io.output import OutputFormat, format_batch_results, format_result

# Console that writes to stderr (so progress doesn't mix with data output)
stderr_console = Console(stderr=True)


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
        if task.total:
            progress = task.completed / task.total
            color_idx = int(progress * len(self.RAINBOW_COLORS) * 3) % len(
                self.RAINBOW_COLORS
            )
        else:
            color_idx = self._color_index
            self._color_index = (self._color_index + 1) % len(self.RAINBOW_COLORS)

        self.complete_style = Style(color=self.RAINBOW_COLORS[color_idx])
        self.finished_style = Style(color="#9400D3")
        return super().render(task)


def create_rainbow_progress() -> Progress:
    """Create a rainbow-colored progress bar that writes to stderr."""
    return Progress(
        SpinnerColumn(style="bold magenta"),
        TextColumn("[bold blue]{task.description}"),
        RainbowBarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=stderr_console,
    )


def get_worker_count(requested: int, num_tasks: int) -> int:
    """Determine the optimal number of workers."""
    if requested == 1 or num_tasks < 10:
        return 1
    cpu_count = os.cpu_count() or 4
    if requested > 0:
        return min(requested, cpu_count)
    return max(1, min(cpu_count - 1, num_tasks))


def run_parallel_batch(
    tasks: list[BatchTask],
    num_workers: int,
    description: str,
) -> tuple[list[WorkerResult], list[str]]:
    """Run batch processing with ProcessPoolExecutor."""
    results: list[WorkerResult] = []
    warnings: list[str] = []

    if num_workers == 1:
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


def find_alignment_files(input_dir: Path) -> list[Path]:
    """Auto-detect alignment files in a directory.

    Tries common FASTA extensions in order: *.fa, *.fasta, *.fna
    Returns files matching the first pattern that finds results.
    """
    for pattern in ["*.fa", "*.fasta", "*.fna"]:
        files = sorted(input_dir.glob(pattern))
        if files:
            return files
    return []


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        typer.echo(f"mkado {__version__}")
        raise typer.Exit()


app = typer.Typer(
    name="mkado",
    help="mkado: McDonald-Kreitman test toolkit.\n\n"
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
    """mkado: McDonald-Kreitman test toolkit."""
    pass


@app.command()
def test(
    fasta: Annotated[
        Path,
        typer.Argument(help="FASTA file (combined alignment, or ingroup sequences)"),
    ],
    outgroup_file: Annotated[
        Optional[Path],
        typer.Argument(help="Outgroup FASTA file (only for separate files mode)"),
    ] = None,
    # === Sequence filtering (combined file mode) ===
    ingroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--ingroup-match", "-i",
            help="Ingroup sequence name pattern (enables combined file mode)",
        ),
    ] = None,
    outgroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--outgroup-match", "-o",
            help="Outgroup sequence name pattern (required with -i)",
        ),
    ] = None,
    polarize_match: Annotated[
        Optional[str],
        typer.Option(
            "--polarize-match",
            help="Second outgroup pattern for polarized test (combined mode)",
        ),
    ] = None,
    # === Polarization (separate files mode) ===
    polarize_file: Annotated[
        Optional[Path],
        typer.Option(
            "--polarize", "-p",
            help="Second outgroup file for polarized test (separate files mode)",
        ),
    ] = None,
    # === Analysis type ===
    use_asymptotic: Annotated[
        bool,
        typer.Option(
            "--asymptotic", "-a",
            help="Use asymptotic MK test (accounts for slightly deleterious mutations)",
        ),
    ] = False,
    # === Asymptotic options ===
    bins: Annotated[
        int,
        typer.Option("--bins", "-b", help="Frequency bins for asymptotic test"),
    ] = 10,
    bootstrap: Annotated[
        int,
        typer.Option("--bootstrap", help="Bootstrap replicates for CI (asymptotic)"),
    ] = 100,
    # === Common options ===
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format: pretty, tsv, json"),
    ] = "pretty",
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1-3)"),
    ] = 1,
    pool_polymorphisms: Annotated[
        bool,
        typer.Option("--pool-polymorphisms", help="Pool polymorphisms from both populations"),
    ] = False,
    min_freq: Annotated[
        float,
        typer.Option("--min-freq", help="Min derived allele frequency (0.0-1.0)"),
    ] = 0.0,
) -> None:
    """Run McDonald-Kreitman test on a single alignment.

    TWO MODES OF OPERATION:

    1. COMBINED FILE MODE (recommended):
       Use -i and -o to filter sequences by name pattern from a single alignment.

       mkado test alignment.fa -i "speciesA" -o "speciesB"
       mkado test alignment.fa -i "gamb" -o "002019" --asymptotic

    2. SEPARATE FILES MODE:
       Provide two FASTA files (ingroup and outgroup).

       mkado test ingroup.fa outgroup.fa
       mkado test ingroup.fa outgroup.fa -p outgroup2.fa  # polarized

    ANALYSIS TYPES:

    - Standard MK test (default): Classic McDonald-Kreitman test
    - Asymptotic MK test (-a): Accounts for slightly deleterious mutations
      by examining alpha across the frequency spectrum (Messer & Petrov 2013)

    EXAMPLES:

        mkado test alignment.fa -i "dmel" -o "dsim"
        mkado test alignment.fa -i "dmel" -o "dsim" -a -b 20
        mkado test alignment.fa -i "dmel" -o "dsim" --polarize-match "dyak"
        mkado test ingroup.fa outgroup.fa
        mkado test ingroup.fa outgroup.fa -a
    """
    from mkado.core.sequences import SequenceSet

    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'.", err=True)
        raise typer.Exit(1)

    fmt = OutputFormat(output_format)

    # Determine mode
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined file mode
        if not ingroup_match or not outgroup_match:
            typer.echo(
                "Error: Combined mode requires both -i/--ingroup-match and -o/--outgroup-match",
                err=True,
            )
            raise typer.Exit(1)

        if outgroup_file is not None:
            typer.echo(
                "Error: Don't provide outgroup file when using -i/-o (combined mode)",
                err=True,
            )
            raise typer.Exit(1)

        if polarize_file is not None:
            typer.echo(
                "Error: Use --polarize-match instead of -p in combined mode",
                err=True,
            )
            raise typer.Exit(1)

        # Load and filter sequences
        all_seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)
        ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
        outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

        if len(ingroup_seqs) == 0:
            typer.echo(f"Error: No sequences match ingroup pattern '{ingroup_match}'", err=True)
            raise typer.Exit(1)
        if len(outgroup_seqs) == 0:
            typer.echo(f"Error: No sequences match outgroup pattern '{outgroup_match}'", err=True)
            raise typer.Exit(1)

        typer.echo(
            f"Found {len(ingroup_seqs)} ingroup, {len(outgroup_seqs)} outgroup sequences",
            err=True,
        )

        # Run appropriate test
        if use_asymptotic:
            if polarize_match:
                typer.echo("Error: Polarized asymptotic test not supported", err=True)
                raise typer.Exit(1)
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
                typer.echo(f"Error: No sequences match polarize pattern '{polarize_match}'", err=True)
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
        # Separate files mode
        if outgroup_file is None:
            typer.echo(
                "Error: Provide outgroup file, or use -i/-o for combined file mode",
                err=True,
            )
            raise typer.Exit(1)

        if polarize_match is not None:
            typer.echo(
                "Error: Use -p/--polarize instead of --polarize-match in separate files mode",
                err=True,
            )
            raise typer.Exit(1)

        # Run appropriate test
        if use_asymptotic:
            if polarize_file:
                typer.echo("Error: Polarized asymptotic test not supported", err=True)
                raise typer.Exit(1)
            result = asymptotic_mk_test(
                ingroup=fasta,
                outgroup=outgroup_file,
                reading_frame=reading_frame,
                num_bins=bins,
                bootstrap_replicates=bootstrap,
                pool_polymorphisms=pool_polymorphisms,
            )
        elif polarize_file:
            result = polarized_mk_test(
                ingroup=fasta,
                outgroup1=outgroup_file,
                outgroup2=polarize_file,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
            )
        else:
            result = mk_test(
                ingroup=fasta,
                outgroup=outgroup_file,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
            )

    typer.echo(format_result(result, fmt))


@app.command()
def batch(
    input_dir: Annotated[
        Path,
        typer.Argument(help="Directory containing FASTA alignment files"),
    ],
    # === Sequence filtering (combined file mode) ===
    ingroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--ingroup-match", "-i",
            help="Ingroup sequence name pattern (enables combined file mode)",
        ),
    ] = None,
    outgroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--outgroup-match", "-o",
            help="Outgroup sequence name pattern (required with -i)",
        ),
    ] = None,
    polarize_match: Annotated[
        Optional[str],
        typer.Option(
            "--polarize-match",
            help="Second outgroup pattern for polarized test",
        ),
    ] = None,
    # === File patterns ===
    file_pattern: Annotated[
        Optional[str],
        typer.Option(
            "--pattern",
            help="File glob pattern (default: auto-detect *.fa, *.fasta, *.fna)",
        ),
    ] = None,
    # === Separate files mode options ===
    ingroup_pattern: Annotated[
        str,
        typer.Option(help="Ingroup file pattern (separate files mode)"),
    ] = "*_ingroup.fa",
    outgroup_pattern: Annotated[
        str,
        typer.Option(help="Outgroup file pattern (separate files mode)"),
    ] = "*_outgroup.fa",
    polarize_pattern: Annotated[
        Optional[str],
        typer.Option(
            "--polarize-pattern",
            help="Second outgroup file pattern (separate files polarized mode)",
        ),
    ] = None,
    # === Analysis type ===
    use_asymptotic: Annotated[
        bool,
        typer.Option(
            "--asymptotic", "-a",
            help="Use asymptotic MK test",
        ),
    ] = False,
    aggregate: Annotated[
        bool,
        typer.Option(
            "--aggregate/--per-gene",
            help="Aggregate across genes (default) or per-gene results (asymptotic)",
        ),
    ] = True,
    # === Asymptotic options ===
    bins: Annotated[
        int,
        typer.Option("--bins", "-b", help="Frequency bins (asymptotic)"),
    ] = 10,
    bootstrap: Annotated[
        int,
        typer.Option("--bootstrap", help="Bootstrap replicates (asymptotic)"),
    ] = 100,
    freq_cutoffs: Annotated[
        str,
        typer.Option("--freq-cutoffs", help="Frequency range 'low,high' (asymptotic)"),
    ] = "0.1,0.9",
    # === Common options ===
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format: pretty, tsv, json"),
    ] = "tsv",
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1-3)"),
    ] = 1,
    pool_polymorphisms: Annotated[
        bool,
        typer.Option("--pool-polymorphisms", help="Pool polymorphisms from both populations"),
    ] = False,
    min_freq: Annotated[
        float,
        typer.Option("--min-freq", help="Min derived allele frequency"),
    ] = 0.0,
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", min=0, help="Parallel workers (0=auto, 1=sequential)"),
    ] = 0,
) -> None:
    """Run MK test on multiple alignment files.

    TWO MODES OF OPERATION:

    1. COMBINED FILE MODE (recommended):
       Each file contains sequences from multiple species.
       Use -i and -o to filter by sequence name pattern.
       Mode is auto-detected when -i is provided.

       mkado batch alignments/ -i "speciesA" -o "speciesB"

    2. SEPARATE FILES MODE:
       Pairs of ingroup/outgroup files (e.g., gene1_ingroup.fa, gene1_outgroup.fa).
       Used when -i is NOT provided.

       mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

    FILE DETECTION:
       Files are auto-detected by trying *.fa, *.fasta, *.fna in order.
       Override with --pattern "*.your_extension"

    ANALYSIS TYPES:

    - Standard MK test (default): Per-gene Dn, Ds, Pn, Ps
    - Asymptotic MK (-a): Genome-wide alpha with frequency correction
      - --aggregate (default): Pool data across genes, fit single curve
      - --per-gene: Fit curve for each gene separately

    EXAMPLES:

        mkado batch alignments/ -i "dmel" -o "dsim"
        mkado batch alignments/ -i "dmel" -o "dsim" -a
        mkado batch alignments/ -i "dmel" -o "dsim" -a --per-gene
        mkado batch alignments/ -i "dmel" -o "dsim" -w 8
        mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"
    """
    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'.", err=True)
        raise typer.Exit(1)

    fmt = OutputFormat(output_format)

    # Parse frequency cutoffs
    try:
        cutoff_parts = freq_cutoffs.split(",")
        frequency_cutoffs = (float(cutoff_parts[0]), float(cutoff_parts[1]))
    except (ValueError, IndexError):
        typer.echo(f"Error: Invalid frequency cutoffs '{freq_cutoffs}'. Use 'low,high'", err=True)
        raise typer.Exit(1)

    # Auto-detect mode based on -i flag
    combined_mode = ingroup_match is not None

    if combined_mode:
        # === COMBINED FILE MODE ===
        if not outgroup_match:
            typer.echo("Error: -o/--outgroup-match required with -i/--ingroup-match", err=True)
            raise typer.Exit(1)

        if use_asymptotic and polarize_match:
            typer.echo("Error: --asymptotic and --polarize-match are mutually exclusive", err=True)
            raise typer.Exit(1)

        # Find alignment files
        if file_pattern:
            alignment_files = sorted(input_dir.glob(file_pattern))
        else:
            alignment_files = find_alignment_files(input_dir)

        if not alignment_files:
            pattern_msg = file_pattern or "*.fa, *.fasta, *.fna"
            typer.echo(f"No files found in {input_dir} (tried: {pattern_msg})", err=True)
            raise typer.Exit(1)

        typer.echo(f"Found {len(alignment_files)} alignment files", err=True)

        # Determine workers
        num_workers = get_worker_count(workers, len(alignment_files))

        # Build tasks
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

        # Aggregated asymptotic mode
        if use_asymptotic and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

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
            tasks, num_workers, "Processing alignments"
        )

        for warning in warnings:
            typer.echo(warning, err=True)

        results = [(r.gene_id, r.result) for r in worker_results]

    else:
        # === SEPARATE FILES MODE ===
        if use_asymptotic and polarize_pattern:
            typer.echo("Error: --asymptotic and --polarize-pattern are mutually exclusive", err=True)
            raise typer.Exit(1)

        ingroup_files = sorted(input_dir.glob(ingroup_pattern))

        if not ingroup_files:
            typer.echo(f"No files matching '{ingroup_pattern}' in {input_dir}", err=True)
            raise typer.Exit(1)

        typer.echo(f"Found {len(ingroup_files)} ingroup files", err=True)

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
                matches = [f for f in potential_matches if f.stem.startswith(prefix)]
                if matches:
                    return matches[0]
                return None
            return outgroup_file_path

        def find_outgroup2_file(ingroup_file: Path) -> Path | None:
            if not polarize_pattern:
                return None
            base_name = ingroup_file.stem
            potential_matches = list(input_dir.glob(polarize_pattern))
            prefix = base_name.split("_")[0]
            matches = [f for f in potential_matches if f.stem.startswith(prefix)]
            return matches[0] if matches else None

        num_workers = get_worker_count(workers, len(ingroup_files))

        tasks = []
        pre_warnings = []
        for ingroup_file in ingroup_files:
            outgroup_file = find_outgroup_file(ingroup_file)
            if outgroup_file is None:
                pre_warnings.append(f"Warning: No outgroup for {ingroup_file.name}")
                continue

            outgroup2_file: Path | None = None
            if polarize_pattern:
                outgroup2_file = find_outgroup2_file(ingroup_file)
                if outgroup2_file is None:
                    pre_warnings.append(f"Warning: No second outgroup for {ingroup_file.name}")
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

        for warning in pre_warnings:
            typer.echo(warning, err=True)

        if not tasks:
            typer.echo("No valid file pairs found", err=True)
            raise typer.Exit(1)

        # Aggregated asymptotic mode
        if use_asymptotic and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

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

        for warning in warnings:
            typer.echo(warning, err=True)

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
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1-3)"),
    ] = 1,
) -> None:
    """Display information about a FASTA file."""
    from mkado.core.sequences import SequenceSet

    seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)

    typer.echo(f"File: {fasta.name}")
    typer.echo(f"Sequences: {len(seqs)}")

    if seqs.sequences:
        typer.echo(f"Alignment length: {seqs.alignment_length} bp")
        typer.echo(f"Codons: {seqs.num_codons}")
        typer.echo(f"Reading frame: {reading_frame}")

        poly_sites = seqs.polymorphic_codons()
        typer.echo(f"Polymorphic codons: {len(poly_sites)}")

        typer.echo("\nSequences:")
        for seq in seqs.sequences:
            typer.echo(f"  {seq.name} ({len(seq)} bp)")


if __name__ == "__main__":
    app()
