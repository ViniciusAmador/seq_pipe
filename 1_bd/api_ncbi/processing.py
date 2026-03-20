import csv
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


def load_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        return list(reader)


def write_csv_rows(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def normalize(value: str) -> str:
    return (value or "").strip()


def normalize_source(source: str) -> str:
    value = normalize(source).casefold()
    if value == "genbank":
        return "genbank"
    if value == "refseq":
        return "refseq"
    return "other"


def filter_rows_with_sequencing_tech(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    return [row for row in rows if normalize(row.get("sequencing_tech", "")) != ""]


def deduplicate_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    has_refseq: Dict[str, bool] = {}
    for row in rows:
        biosample = normalize(row.get("biosample", ""))
        source = normalize_source(row.get("source", ""))
        if biosample not in has_refseq:
            has_refseq[biosample] = False
        if source == "refseq":
            has_refseq[biosample] = True

    treated: List[Dict[str, str]] = []
    for row in rows:
        biosample = normalize(row.get("biosample", ""))
        source = normalize_source(row.get("source", ""))
        if has_refseq.get(biosample, False) and source != "refseq":
            continue
        treated.append(row)
    return treated


def apply_dedup_pipeline(rows: Sequence[Dict[str, str]]) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    filtered = filter_rows_with_sequencing_tech(rows)
    treated = deduplicate_rows(filtered)

    def source_counter(data_rows: Sequence[Dict[str, str]]) -> Dict[str, int]:
        counter = Counter(normalize_source(row.get("source", "")) for row in data_rows)
        return {
            "genbank": counter.get("genbank", 0),
            "refseq": counter.get("refseq", 0),
            "other": counter.get("other", 0),
        }

    stats = {
        "rows_before": len(rows),
        "rows_removed_missing_sequencing_tech": len(rows) - len(filtered),
        "rows_after_sequencing_filter": len(filtered),
        "rows_removed_non_refseq_when_biosample_has_refseq": len(filtered) - len(treated),
        "rows_after_dedup": len(treated),
        "source_before_genbank": source_counter(rows)["genbank"],
        "source_before_refseq": source_counter(rows)["refseq"],
        "source_after_seq_filter_genbank": source_counter(filtered)["genbank"],
        "source_after_seq_filter_refseq": source_counter(filtered)["refseq"],
        "source_after_dedup_genbank": source_counter(treated)["genbank"],
        "source_after_dedup_refseq": source_counter(treated)["refseq"],
    }
    return treated, stats
