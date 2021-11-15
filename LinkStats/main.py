#!/usr/bin/env python

import sys
import warnings
from concurrent.futures import ThreadPoolExecutor as TPE
from dataclasses import dataclass
from enum import Enum, auto
from itertools import chain, groupby, tee
from pathlib import Path
from typing import List

import click as ck
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pysam
from pysam import AlignmentFile as AF
from tqdm import tqdm
from tqdm.contrib.concurrent import thread_map as tqdm_map

from _LinkStats_C import _LinkStats

pysam.set_verbosity(0)


warnings.filterwarnings("ignore", message="divide by zero")
warnings.filterwarnings("ignore", message="invalid value encountered in double_scalars")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")


mpl.use("agg")
mpl.rc("font", **{"family": "sans", "weight": "normal", "size": 14})
mpl.rcParams["agg.path.chunksize"] = 10000
plt.style.use("ggplot")
plt.rcParams["figure.figsize"] = (10.0, 7.0)


import seaborn as sb

sb.set(style="darkgrid", color_codes=True)


class CallBack:
    class Types(Enum):
        AF = auto()
        MD = auto()
        HD = auto()
        CP = auto()
        PP = auto()

    def __init__(self, cb_type):
        self.cb_type = cb_type

    @property
    def is_AF(self):
        return self.cb_type == self.Types.AF

    @property
    def is_MD(self):
        return self.cb_type == self.Types.MD

    @property
    def is_HD(self):
        return self.cb_type == self.Types.HD

    @property
    def is_CP(self):
        return self.cb_type == self.Types.CP

    @property
    def is_PP(self):
        return self.cb_type == self.Types.PP


class AlignmentFile(CallBack):
    def __init__(self, file, ref, use_mi_tags, cluster_threshold, name=None):
        super().__init__(self.Types.AF)

        self.file = file
        self.ref = ref
        self.use_mi_tags = use_mi_tags
        self.cluster_threshold = cluster_threshold
        self.name = name

    @property
    def stem_name(self):
        return self.file.stem if str(self.file) != "-" else "<stdin>"


class Data(CallBack):
    def __init__(self, file, cb_type):
        super().__init__(cb_type)

        self.file = file


class MoleculeData(Data):
    def __init__(self, file):
        super().__init__(file, self.Types.MD)


class HistData(Data):
    def __init__(self, file):
        super().__init__(file, self.Types.HD)


class Prefix(CallBack):
    def __init__(self, prefix, cb_type):
        super().__init__(cb_type)

        self.prefix = prefix


class CSVPrefix(Prefix):
    def __init__(self, prefix, save_summ, save_mol, save_hist):
        super().__init__(prefix, self.Types.CP)
        self.save_summ = save_summ
        self.save_mol = save_mol
        self.save_hist = save_hist

    @property
    def any_set(self):
        return self.save_summ or self.save_mol or self.save_hist


class PlotPrefix(Prefix):
    def __init__(self, prefix):
        super().__init__(prefix, self.Types.PP)


def ConcatDF(it):
    def iter_then_empty():
        yield from it
        yield pd.DataFrame()

    return pd.concat(iter_then_empty(), ignore_index=True)


def GetAllStats(alignment_files, molecular_data, min_reads, threads):
    alignment_files = tuple(alignment_files)

    def GetAllStatsFromAFs():
        @dataclass(frozen=True, eq=True)
        class Alignment:
            haveMI: bool
            mapping_quality: int
            reference_start: int
            reference_end: int
            query_length: int
            mi: int

        class Molecule:
            def __init__(self, alignment, bx, reference_name):
                self.min = np.inf
                self.max = 0

                self.n_reads = 0
                self.total_read_len = 0
                self.total_mapping_quality = 0

                self.mi = alignment.mi if alignment.haveMI else None
                self.bx = bx
                self.ref = reference_name

            def update(self, alignment):
                self.min = min(self.min, alignment.reference_start)
                self.max = max(self.max, alignment.reference_end)

                self.n_reads += 1
                self.total_read_len += alignment.query_length
                self.total_mapping_quality += alignment.mapping_quality

            def __len__(self):
                return max(0, self.max - self.min)

            @property
            def mean_mapq(self):
                return self.total_mapping_quality / self.n_reads

            @property
            def mean_read_depth(self):
                return self.total_read_len / len(self)

        @dataclass(frozen=True, eq=True)
        class BasicStats:
            insert_sizes: List[int]
            total_read_length: int
            total_alignments: int
            total_dup: int
            total_qcf: int
            total_unm: int
            total_nomi: int
            total_nobx: int
            total_zeromq: int

        def ReadSAM(alignment_file, threads):
            @dataclass
            class _LinkStats_Params:
                log: int
                num_threads: int
                sam_file_name: str
                fasta_file_name: str
                override_name: str
                fallback_name: str
                use_mi: bool

            genome_length, ref_names, basic_stats, molecule_data = _LinkStats(
                _LinkStats_Params(
                    log=sys.stderr,
                    num_threads=threads,
                    sam_file_name=str(alignment_file.file),
                    fasta_file_name=alignment_file.ref,
                    override_name=alignment_file.name,
                    fallback_name=alignment_file.stem_name,
                    use_mi=alignment_file.use_mi_tags,
                )
            )

            return (
                genome_length,
                {
                    name: BasicStats(list(stats[0]), *stats[1:])
                    for name, stats in basic_stats
                },
                {
                    name: {
                        ref_names[tid]: {
                            bxmi: tuple(Alignment(*aln) for aln in alns)
                            for bxmi, alns in t2
                        }
                        for tid, t2 in t1
                    }
                    for name, t1 in molecule_data
                },
            )

        def GetStats(alignment_file, threads):
            genome_length, basic_stats, molecule_data = ReadSAM(alignment_file, threads)

            def ClusterAlignments(alignments, bx, reference_name):
                def pairwise(it):
                    a, b = tee(it)
                    next(b, None)
                    return zip(a, b)

                g = nx.Graph()
                g.add_nodes_from(alignments)
                g.add_edges_from(
                    (a, b)
                    for a, b in pairwise(
                        sorted(alignments, key=lambda a: a.reference_start)
                    )
                    if (b.reference_start - a.reference_end)
                    < alignment_file.cluster_threshold
                )
                for cc in nx.connected_components(g):
                    cc = iter(cc)
                    c = next(cc)
                    m = Molecule(c, bx, reference_name)
                    for c in chain((c,), cc):
                        m.update(c)
                    yield m

            molecule_data = pd.DataFrame(
                (
                    (
                        name,
                        len(molecule),
                        molecule.n_reads,
                        molecule.mi,
                        molecule.bx,
                        molecule.ref,
                        molecule.mean_read_depth,
                        molecule.mean_mapq,
                    )
                    for name, a in molecule_data.items()
                    for reference_name, b in a.items()
                    for (bx, _), c in b.items()
                    for molecule in ClusterAlignments(c, bx, reference_name)
                ),
                columns=(
                    "Sample Name",
                    "Molecule Length",
                    "No. Reads",
                    "MI",
                    "BX",
                    "Reference",
                    "Mean Read Depth",
                    "Mean MapQ",
                ),
            )

            def n_stats(data, ns):
                d = np.cumsum((0,) + tuple(np.sort(data)[::-1]))
                return tuple(
                    d[x] - d[x - 1]
                    for n in ns
                    for x in (np.where(d >= (d[-1] * n / 100))[0].min(),)
                )

            return (
                molecule_data,
                pd.DataFrame(
                    (
                        (
                            (
                                name,
                                genome_length,
                                bs.total_alignments,
                                bs.total_dup / bs.total_alignments,
                                bs.total_qcf / bs.total_alignments,
                                bs.total_unm / bs.total_alignments,
                                (bs.total_alignments - bs.total_unm)
                                / bs.total_alignments,
                                bs.total_nobx / bs.total_alignments,
                                bs.total_nomi / bs.total_alignments,
                                bs.total_zeromq / bs.total_alignments,
                            )
                            + n_stats(data["No. Reads"], (50, 90))
                            + (
                                (data["No. Reads"] ** 2).sum()
                                / data["No. Reads"].sum(),
                            )
                            + tuple(
                                chain.from_iterable(
                                    (
                                        (
                                            (
                                                q.shape[0],
                                                q.mean(),
                                                d.mean(),
                                            )
                                            + n_stats(
                                                d,
                                                (50, 90),
                                            )
                                            + ((d ** 2).sum() / d.sum(),)
                                        )
                                        for m in min_reads
                                        for s in (data[data["No. Reads"] >= m],)
                                        for d in (s["Molecule Length"],)
                                        for q in (s["Mean MapQ"],)
                                    )
                                )
                            )
                            + (
                                np.median(bs.insert_sizes),
                                bs.total_read_length / genome_length,
                            )
                            + tuple(
                                chain.from_iterable(
                                    (
                                        (
                                            s["Mean Read Depth"].mean(),
                                            s["Molecule Length"].sum() / genome_length,
                                        )
                                        for m in min_reads
                                        for s in (data[data["No. Reads"] >= m],)
                                    )
                                )
                            )
                        )
                        for name in molecule_data["Sample Name"].unique()
                        for bs in (basic_stats[name],)
                        for data in (
                            molecule_data[molecule_data["Sample Name"] == name],
                        )
                    ),
                    columns=(
                        "Sample Name",
                        "Genome Length",
                        "Total Alignments",
                        "Duplicates",
                        "QCFail",
                        "Unmapped",
                        "Mapped",
                        "No BX",
                        "No MI",
                        "Zero MapQ",
                        "N50 Reads Per Molecule",
                        "N90 Reads Per Molecule",
                        "auN Reads Per Molecule",
                    )
                    + tuple(
                        chain.from_iterable(
                            (
                                (
                                    f"No. Molecules (No. Reads >= {m})",
                                    f"Mean Read MapQ Per Molecule (No. Reads >= {m})",
                                    f"Mean Molecule Length (No. Reads >= {m})",
                                    f"N50 Molecule Length (No. Reads >= {m})",
                                    f"N90 Molecule Length (No. Reads >= {m})",
                                    f"auN Molecule Length (No. Reads >= {m})",
                                )
                                for m in min_reads
                            )
                        )
                    )
                    + ("Median Insert Size", "Mean Short Read Depth")
                    + tuple(
                        chain.from_iterable(
                            (
                                (
                                    f"Mean Short Read Depth Per Molecule (No. Reads >= {m})",
                                    f"Molecule Read Depth (No. Reads >= {m})",
                                )
                                for m in min_reads
                            )
                        )
                    ),
                ),
            )

        max_workers = min(threads, len(alignment_files))
        with TPE(max_workers=max_workers) as exe:
            return exe.map(
                lambda af: GetStats(af, threads=max(threads // max_workers, 1)),
                alignment_files,
            )

    def GetAllStatsFromCSVs():
        mol_files = tuple(mol.file for mol in molecular_data)

        return (
            iter(
                tqdm_map(
                    pd.read_csv,
                    mol_files,
                    max_workers=threads,
                    desc="Read CSVs (molecular data)",
                    unit=" CSV files",
                    unit_scale=True,
                )
            )
            if len(mol_files) > 0
            else ()
        )

    summary_dfs = []

    def yield_all():
        for df, summ_df in GetAllStatsFromAFs():
            summary_dfs.append(summ_df)
            yield df
        yield from GetAllStatsFromCSVs()

    return ConcatDF(yield_all()), ConcatDF(summary_dfs)


def GetAllMolLenHists(df, hist_data, min_reads, threads):
    def GetMolLenHist(args):
        MAX_BINS = 1024

        sample_name, min_reads = args
        data = df[(df["Sample Name"] == sample_name) & (df["No. Reads"] >= min_reads)][
            "Molecule Length"
        ]
        prob, length = np.histogram(
            data,
            bins=np.interp(
                np.linspace(
                    0,
                    len(data),
                    np.clip(
                        len(np.histogram_bin_edges(data, bins="auto")) - 1, 1, MAX_BINS
                    )
                    + 1,
                ),
                np.arange(len(data)),
                np.sort(data),
            ),
            density=True,
        )

        select = ~np.isnan(prob)

        return pd.DataFrame(
            {
                "PDF": prob[select],
                "CDF": np.cumsum(prob[select]) / prob[select].sum(),
                "Molecule Length": ((length[:-1] + length[1:]) / 2)[select],
                "Sample Name": sample_name,
                "Min No. Reads": str(min_reads),
            }
        )

    def yield_all():
        if df.shape[0] > 0:
            yield from iter(
                tqdm_map(
                    GetMolLenHist,
                    tuple(
                        (name, n)
                        for name in df["Sample Name"].unique()
                        for n in min_reads
                    ),
                    max_workers=threads,
                    desc="Generate Histogram Data",
                    unit=" Data-Sets",
                    unit_scale=True,
                )
            )

        hist_files = tuple(hist.file for hist in hist_data)

        yield from (
            iter(
                tqdm_map(
                    pd.read_csv,
                    hist_files,
                    max_workers=threads,
                    desc="Read CSVs (histogram data)",
                    unit=" CSV files",
                    unit_scale=True,
                )
            )
            if len(hist_files) > 0
            else ()
        )

    return ConcatDF(yield_all())


def documenter(docstring):
    def inner_documenter(f):
        f.__doc__ = docstring
        return f

    return inner_documenter


@ck.group(chain=True)
@ck.option(
    "-t",
    "--threads",
    type=ck.IntRange(1, None, clamp=True),
    default=4,
    help="Number of threads to use, Default=4",
)
@ck.option(
    "-m",
    "--min_reads",
    type=ck.IntRange(1, None, clamp=True),
    multiple=True,
    default=(1, 3, 5, 10),
    help="Minimum reads per molecule for analysis, multiple values possible, Default=(1, 3, 5, 10)",
)
@ck.version_option()
@documenter(
    """
Collect summary and molecular data from aligned short linked-reads in SAM format.

\b
Combine multiple data sources into summary plots.

\b
Usage Example, read SAM/BAM/CRAM from <stdin> and save the summary and molecule data in csv format. Analyse molecules grouped by 5 and 10 minimum reads per molecule.
-------------
...<sam/bam/cram> | LinkStats -t 16 -m 5 -m 10 sam-data - save-csvs results/csvs/

\b
Usage Example, combine histogram data from multiple sources into summary plots.
-------------
LinkStats -t 16 hist-data results/dataset_1_molecular_length_histograms.csv.bz2 hist-data results/dataset_2_molecular_length_histograms.csv.bz2 hist-data results/dataset_3_molecular_length_histograms.csv.bz2 save-plots results/plots/
"""
)
def cli(threads, min_reads):
    pass


@cli.command()
@ck.argument("path", type=ck.Path(readable=True, path_type=Path))
@ck.option("-r", "--reference", type=ck.Path(exists=True), help="FASTA reference")
@ck.option("-n", "--name", type=str, help="Sample name")
@ck.option("--mi/--no-mi", default=False, help="Group by MI tags as well as BX")
@ck.option(
    "-t",
    "--threshold",
    type=int,
    default=50000,
    help="Maximum alignment separation per molecule",
)
@documenter(
    """
Read SAM/BAM/CRAM data from PATH.

\b
Creates summary and molecular data-sets for each sample-name (SM tag).

\b
Alignments must have MI (molecular id) and BX (barcode) tags.

\b
Option -n sets a same-name for all data, taking preference over SM tags.

\b
Option -r sets a FASTA reference for CRAM decoding.
"""
)
def sam_data(path, mi, threshold, reference=None, name=None):
    return AlignmentFile(
        file=path, ref=reference, name=name, use_mi_tags=mi, cluster_threshold=threshold
    )


@cli.command()
@ck.argument("file", type=ck.Path(readable=True, path_type=Path))
@documenter(
    """
Read in molecular data from a CSV FILE.

\b
Use to re-calculate histogram data.
"""
)
def molecule_data(file):
    return MoleculeData(file=file)


@cli.command()
@ck.argument("file", type=ck.Path(readable=True, path_type=Path))
@documenter(
    """
Read in molecule length histogram data from a CSV FILE.

\b
Use to re-generate or create combined plots.
"""
)
def hist_data(file):
    return HistData(file=file)


@cli.command()
@ck.argument("prefix", type=ck.Path(readable=True, path_type=Path))
@ck.option("--summ/--no-summ", default=True, help="Save summary data table")
@ck.option("--mol/--no-mol", default=False, help="Save molecule data table")
@ck.option("--hist/--no-hist", default=False, help="Save histogram data table")
@documenter(
    """
Saves summary, molecule or histogram data to CSV files at PREFIX_

\b
By default, only summary data is saved
"""
)
def save_csvs(prefix, summ, mol, hist):
    return CSVPrefix(prefix, save_summ=summ, save_mol=mol, save_hist=hist)


@cli.command()
@ck.argument("prefix", type=ck.Path(readable=True, path_type=Path))
@documenter(
    """
Generates plots from any histogram data and saves them at PREFIX_
"""
)
def save_plots(prefix):
    return PlotPrefix(prefix)


@cli.result_callback()
def run(callbacks, threads, min_reads):
    print("Starting...\n", file=sys.stderr)

    csv = tuple(cp for cp in callbacks if cp.is_CP)
    if len(csv) == 0:
        csv = None
    else:
        if len(csv) > 1:
            warnings.warn(
                f"More than one CSV prefix specified, using last one: {csv[-1].prefix}"
            )
        csv = csv[-1]

    if csv and not csv.any_set:
        raise ck.ClickException("CSV prefix specified, but no data set to be saved")

    plot = tuple(pp.prefix for pp in callbacks if pp.is_PP)
    if len(plot) == 0:
        plot = None
    else:
        if len(plot) > 1:
            warnings.warn(
                f"More than one Plot prefix specified, using last one: {plot[-1]}"
            )
        plot = plot[-1]

    if not (csv or plot):
        raise ck.ClickException("Neither CSV nor Plot prefix specified, nothing to do")

    all_data, summary_data = GetAllStats(
        (af for af in callbacks if af.is_AF),
        (md for md in callbacks if md.is_MD),
        min_reads,
        threads,
    )

    if summary_data.shape[0] > 0:
        print("\n", summary_data, "\n", sep="", file=sys.stderr)

    hist_data = GetAllMolLenHists(
        all_data, (hd for hd in callbacks if hd.is_HD), min_reads, threads
    )

    def base_get_path(f, prefix):
        return prefix / f if prefix.is_dir() else Path(str(prefix) + "_" + f)

    generated = []

    if csv:
        csv.prefix.parent.mkdir(parents=True, exist_ok=True)
        get_path = lambda f: base_get_path(f, csv.prefix)

        def save_csv(args):
            df, name = args
            df.to_csv(get_path(name), index=False)
            return get_path(name)

        generated.append(
            iter(
                tqdm_map(
                    save_csv,
                    (
                        (((summary_data, "summary_data.csv"),) if csv.save_summ else ())
                        + (
                            ((all_data, "molecular_data.csv.bz2"),)
                            if csv.save_mol
                            else ()
                        )
                        + (
                            ((hist_data, "molecular_length_histograms.csv.bz2"),)
                            if csv.save_hist
                            else ()
                        )
                    ),
                    max_workers=threads,
                    desc="Saving CSV data",
                    unit=" Data-Sets",
                )
            )
        )

    if plot:
        plot.parent.mkdir(parents=True, exist_ok=True)
        get_path = lambda f: base_get_path(f, plot)

        def save_plots(col, hue, n, typ):
            name = get_path(f"molecular_length_{typ}s_{n}.png")
            sb.relplot(
                kind="line",
                data=hist_data,
                col=col,
                hue=hue,
                x="Molecule Length",
                y=typ,
            ).set(xscale="log", yscale=("log" if typ == "PDF" else "linear")).savefig(
                name,
                dpi=200,
                bbox_inches="tight",
            )
            return name

        generated.append(
            (
                save_plots(col, hue, i + 1, typ)
                for i, col, hue, typ in tqdm(
                    tuple(
                        (i, col, hue, typ)
                        for i, (col, hue) in enumerate(
                            (col, hue)
                            for colhue in (("Sample Name", "Min No. Reads"),)
                            for col, hue in (colhue, colhue[::-1])
                        )
                        for typ in ("PDF", "CDF")
                    ),
                    desc="Saving Plots",
                    unit=" Plots",
                )
            )
        )

    print(
        "\nGenerated files:",
        *(("\t" + str(n)) for n in chain.from_iterable(generated)),
        "\nDone",
        sep="\n",
        file=sys.stderr,
    )


if __name__ == "__main__":
    cli()
