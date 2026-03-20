import argparse
import csv
import sys
from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from processing import apply_dedup_pipeline, normalize, normalize_source

FIELDS = ["assembly_accession", "bioproject", "biosample"]
REQUIRED_COLUMNS = FIELDS + ["source"]
SOURCE_SET_ORDER = [
    ("bioproject", "genbank"),
    ("bioproject", "refseq"),
    ("biosample", "genbank"),
    ("biosample", "refseq"),
    ("assembly_accession", "genbank"),
    ("assembly_accession", "refseq"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Gera graficos (pre-dedup) com notacao matematica "
            "para presenca de GenBank/RefSeq em ASM, BP e BS."
        )
    )
    parser.add_argument(
        "--input",
        default="outputs/metadata/assemblies.csv",
        help="CSV antes da deduplicacao (padrao: outputs/metadata/assemblies.csv).",
    )
    parser.add_argument(
        "--outdir",
        default="outputs/metadata",
        help="Diretorio de saida dos graficos (padrao: outputs/metadata).",
    )
    return parser.parse_args()


def resolve_path(base_dir: Path, path_str: str) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path
    if path.exists():
        return path.resolve()
    return base_dir / path


def normalize_key(_field: str, value: str) -> str:
    return normalize(value)


def source_set_key(field: str, source: str) -> str:
    return f"{field}|{source}"


def source_set_pretty_name(set_key: str) -> str:
    field, source = set_key.split("|", 1)
    field_symbol = {
        "assembly_accession": "ASM",
        "bioproject": "BP",
        "biosample": "BS",
    }.get(field, field.upper())
    source_symbol = "GB" if source == "genbank" else "RS"
    return f"{field_symbol}_{source_symbol}"


def source_set_short_name(set_key: str) -> str:
    field, source = set_key.split("|", 1)
    field_short = {
        "assembly_accession": "ASM",
        "bioproject": "BP",
        "biosample": "BS",
    }.get(field, field.upper())
    source_short = "GB" if source == "genbank" else "RS"
    return f"{field_short}_{source_short}"


def build_region_counts(
    row_sets: Dict[str, Set[int]], total_rows: int, set_order: Sequence[str]
) -> Dict[Tuple[bool, ...], int]:
    regions: Dict[Tuple[bool, ...], int] = {}
    for flags in product([False, True], repeat=len(set_order)):
        if any(flags):
            regions[flags] = 0

    for idx in range(total_rows):
        flags = tuple(idx in row_sets[name] for name in set_order)
        if any(flags):
            regions[flags] += 1
    return regions


def combo_short_label(flags: Tuple[bool, ...], set_order: Sequence[str]) -> str:
    active = [
        source_set_short_name(name)
        for name, present in zip(set_order, flags)
        if present
    ]
    return " n ".join(active)


def load_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise SystemExit(f"CSV sem cabecalho: {path}")

        missing = [field for field in REQUIRED_COLUMNS if field not in reader.fieldnames]
        if missing:
            raise SystemExit(
                f"CSV {path} sem colunas obrigatorias: {', '.join(missing)}"
            )
        return list(reader)


def simulate_dedup_counts(rows: Sequence[Dict[str, str]]) -> Dict[str, int]:
    _treated, stats = apply_dedup_pipeline(rows)
    return stats


def analyze_rows(rows: Sequence[Dict[str, str]]) -> Dict[str, object]:
    row_sets: Dict[str, Set[int]] = {}
    source_row_sets: Dict[str, Set[int]] = {}
    stats: Dict[str, Dict[str, int]] = {}

    for field in FIELDS:
        key_sources: Dict[str, Set[str]] = defaultdict(set)
        row_keys: List[str] = []

        for row in rows:
            key = normalize_key(field, row.get(field, ""))
            row_keys.append(key)
            if not key:
                continue
            source = normalize_source(row.get("source", ""))
            if source in {"genbank", "refseq"}:
                key_sources[key].add(source)

        shared_keys = {
            key
            for key, sources in key_sources.items()
            if "genbank" in sources and "refseq" in sources
        }
        only_genbank_keys = {
            key for key, sources in key_sources.items() if sources == {"genbank"}
        }
        only_refseq_keys = {
            key for key, sources in key_sources.items() if sources == {"refseq"}
        }

        row_set_shared = {
            idx for idx, key in enumerate(row_keys) if key and key in shared_keys
        }
        row_set_only_genbank = {
            idx for idx, key in enumerate(row_keys) if key and key in only_genbank_keys
        }
        row_set_only_refseq = {
            idx for idx, key in enumerate(row_keys) if key and key in only_refseq_keys
        }

        row_set_source_genbank = {
            idx
            for idx, key in enumerate(row_keys)
            if key and "genbank" in key_sources.get(key, set())
        }
        row_set_source_refseq = {
            idx
            for idx, key in enumerate(row_keys)
            if key and "refseq" in key_sources.get(key, set())
        }

        row_sets[field] = row_set_shared
        source_row_sets[source_set_key(field, "genbank")] = row_set_source_genbank
        source_row_sets[source_set_key(field, "refseq")] = row_set_source_refseq

        stats[field] = {
            "shared_key_count": len(shared_keys),
            "rows_in_shared_keys": len(row_set_shared),
            "only_genbank_key_count": len(only_genbank_keys),
            "only_refseq_key_count": len(only_refseq_keys),
            "rows_in_only_genbank_keys": len(row_set_only_genbank),
            "rows_in_only_refseq_keys": len(row_set_only_refseq),
            "rows_in_genbank_or_both_keys": len(row_set_source_genbank),
            "rows_in_refseq_or_both_keys": len(row_set_source_refseq),
        }

    field_regions = build_region_counts(row_sets, len(rows), FIELDS)
    source_set_order = [source_set_key(field, source) for field, source in SOURCE_SET_ORDER]
    source_regions = build_region_counts(source_row_sets, len(rows), source_set_order)

    rows_with_any_shared = len(set().union(*row_sets.values()))

    return {
        "total_rows": len(rows),
        "row_sets": row_sets,
        "regions": field_regions,
        "source_row_sets": source_row_sets,
        "source_set_order": source_set_order,
        "source_regions": source_regions,
        "rows_with_any_shared": rows_with_any_shared,
        "stats": stats,
        "dedup_counts": simulate_dedup_counts(rows),
    }


def plot_upset(data: Dict[str, object], title: str, out_png: Path) -> None:
    regions: Dict[Tuple[bool, ...], int] = data["source_regions"]  # type: ignore[assignment]
    row_sets: Dict[str, Set[int]] = data["source_row_sets"]  # type: ignore[assignment]
    set_order: List[str] = data["source_set_order"]  # type: ignore[assignment]
    total_rows: int = data["total_rows"]  # type: ignore[assignment]
    dedup_counts: Dict[str, int] = data["dedup_counts"]  # type: ignore[assignment]

    combos = [(flags, count) for flags, count in regions.items() if count > 0]
    combos.sort(key=lambda item: (-item[1], -sum(item[0]), combo_short_label(item[0], set_order)))
    combos = combos[:20]

    n_sets = len(set_order)
    y_pos = list(range(n_sets))
    x_pos = list(range(len(combos)))

    fig = plt.figure(figsize=(16, 9))
    grid = fig.add_gridspec(2, 2, width_ratios=[2.1, 3.9], height_ratios=[3.0, 1.8])

    ax_blank = fig.add_subplot(grid[0, 0])
    ax_blank.axis("off")

    ax_inter = fig.add_subplot(grid[0, 1])
    counts = [count for _, count in combos]
    ax_inter.bar(x_pos, counts, color="#2c7fb8")
    ax_inter.set_ylabel("|C_i|")
    ax_inter.set_xticks([])
    ax_inter.set_title(r"Top $|C_i|$")
    for x, value in zip(x_pos, counts):
        if value > 0:
            ax_inter.text(x, value + 0.4, str(value), ha="center", va="bottom", fontsize=9)

    ax_set = fig.add_subplot(grid[1, 0])
    set_sizes = [len(row_sets[name]) for name in set_order]
    ax_set.barh(y_pos, set_sizes, color="#74a9cf")
    ax_set.set_yticks(y_pos)
    ax_set.set_yticklabels([source_set_pretty_name(name) for name in set_order])
    ax_set.invert_yaxis()
    ax_set.set_xlabel("|S_j|")
    for y, value in zip(y_pos, set_sizes):
        ax_set.text(value + 0.5, y, str(value), va="center", fontsize=9)

    ax_matrix = fig.add_subplot(grid[1, 1])
    for x, (flags, _) in zip(x_pos, combos):
        active = [idx for idx, is_on in enumerate(flags) if is_on]
        for y in y_pos:
            if flags[y]:
                ax_matrix.scatter(x, y, s=80, color="#045a8d")
            else:
                ax_matrix.scatter(x, y, s=30, color="#d0d0d0")
        if len(active) > 1:
            ax_matrix.plot([x, x], [min(active), max(active)], color="#045a8d", linewidth=2)

    ax_matrix.set_yticks(y_pos)
    ax_matrix.set_yticklabels([])
    ax_matrix.invert_yaxis()
    ax_matrix.set_xticks(x_pos)
    ax_matrix.set_xticklabels(
        [combo_short_label(flags, set_order) for flags, _ in combos],
        rotation=90,
        fontsize=8,
    )
    ax_matrix.set_xlabel("Intersecoes")

    title_lines = [
        title,
        rf"$N_{{rows}}={total_rows}$ | $N_{{shared}}={data['rows_with_any_shared']}$",
        (
            r"$N_{seq\_filter}="
            + str(dedup_counts["rows_after_sequencing_filter"])
            + r"$ | $N_{dedup}="
            + str(dedup_counts["rows_after_dedup"])
            + r"$"
        ),
    ]
    fig.suptitle("\n".join(title_lines), fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def write_summary_csv(data: Dict[str, object], out_csv: Path) -> None:
    stats: Dict[str, Dict[str, int]] = data["stats"]  # type: ignore[assignment]
    dedup_counts: Dict[str, int] = data["dedup_counts"]  # type: ignore[assignment]
    total_rows: int = data["total_rows"]  # type: ignore[assignment]
    rows_with_any_shared: int = data["rows_with_any_shared"]  # type: ignore[assignment]

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerow(["total_rows_pre_dedup", total_rows])
        writer.writerow(["rows_with_any_shared_relation", rows_with_any_shared])
        for field in FIELDS:
            for key, value in stats[field].items():
                writer.writerow([f"{field}.{key}", value])
        for key, value in dedup_counts.items():
            writer.writerow([f"dedup.{key}", value])


def main() -> None:
    args = parse_args()
    base_dir = Path(__file__).resolve().parents[1]
    input_path = resolve_path(base_dir, args.input)
    outdir = resolve_path(base_dir, args.outdir)

    rows = load_rows(input_path)
    data = analyze_rows(rows)

    summary_csv = outdir / "dedup_summary.csv"
    venn_png = outdir / "dedup_story_before.png"
    upset_png = outdir / "dedup_upset_before.png"

    write_summary_csv(data, summary_csv)
    plot_upset(data, "Antes da deduplicacao", upset_png)
    plot_upset(data, "Antes da deduplicacao", venn_png)

    print(f"OK: {summary_csv}")
    print(f"OK: {venn_png}")
    print(f"OK: {upset_png}")


if __name__ == "__main__":
    main()
