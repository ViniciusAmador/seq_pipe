import csv
import json
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import yaml
from Bio import SeqIO
from rich.console import Console, Group
from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.table import Table

from settings import command_available, load_local_dotenv, resolve_ncbi_credentials

console = Console()


@dataclass
class Row:
    sequencing_tech: str
    geo_loc_name: str
    country: str
    assembly_accession: str
    assembly_submission_date: str
    assembly_release_date: str
    assembly_seq_release_date: str
    seed_year: str
    organism_name: str
    taxid: str
    bioproject: str
    biosample: str
    source: str
    assembly_level: str
    genome_size: str
    total_ungapped_length: str
    gaps_between_scaffolds: str
    number_of_chromosomes: str
    number_of_organelles: str
    number_of_scaffolds: str
    scaffold_n50: str
    scaffold_l50: str
    number_of_contigs: str
    contig_n50: str
    contig_l50: str
    gc_percent: str
    observed_genome_size: str
    observed_contig_count: str
    observed_contig_n50: str
    observed_contig_l50: str
    observed_gc_percent: str
    observed_n_count: str
    observed_min_contig_len: str
    observed_max_contig_len: str


@dataclass
class ItemStatus:
    accession: str
    source: str
    zip_ready: bool = False
    fasta_ready: bool = False
    report_ready: bool = False
    state: str = "pending"
    note: str = "aguardando"


@dataclass
class DashboardState:
    job_name: str
    taxon_value: str
    total_hits: int
    statuses: List[ItemStatus] = field(default_factory=list)
    current_accession: str = "-"
    completed: int = 0
    failures: int = 0

    def counters_by_source(self) -> Counter:
        counter: Counter = Counter()
        for item in self.statuses:
            if item.state == "done":
                counter[item.source] += 1
        return counter

    def pending_count(self) -> int:
        return sum(1 for item in self.statuses if item.state == "pending")


# -----------------------------
# Helpers
# -----------------------------
def run_cmd(cmd: List[str], *, env: Optional[Dict[str, str]] = None) -> subprocess.CompletedProcess:
    res = subprocess.run(cmd, capture_output=True, env=env)
    if res.returncode != 0:
        raise RuntimeError(
            "Falha ao executar comando:\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{res.stdout.decode('utf-8', errors='replace')}\n"
            f"STDERR:\n{res.stderr.decode('utf-8', errors='replace')}\n"
        )
    return res


def parse_jsonl_stdout(stdout_bytes: bytes) -> List[Dict[str, Any]]:
    if not stdout_bytes:
        return []
    text = stdout_bytes.decode("utf-8", errors="replace")
    lines = [line for line in text.splitlines() if line.strip()]
    return [json.loads(line) for line in lines]


def safe_get(data: Dict[str, Any], path: List[str], default=None):
    current: Any = data
    for key in path:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def infer_source_from_accession(accession: str) -> str:
    if accession.startswith("GCF_"):
        return "RefSeq"
    if accession.startswith("GCA_"):
        return "GenBank"
    return "Unknown"


def ensure_dirs(base: Path) -> Dict[str, Path]:
    fasta_dir = base / "fasta"
    meta_dir = base / "metadata"
    pkgs_dir = base / "packages"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)
    pkgs_dir.mkdir(parents=True, exist_ok=True)
    return {"fasta": fasta_dir, "meta": meta_dir, "pkgs": pkgs_dir}


def is_valid_zip(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    import zipfile

    return zipfile.is_zipfile(path)


def run_cmd_with_retries(
    cmd: List[str],
    *,
    retries: int,
    backoff_seconds: int,
    env: Optional[Dict[str, str]] = None,
) -> subprocess.CompletedProcess:
    retryable = [
        "timeout",
        "timed out",
        "tls handshake",
        "connection",
        "temporary",
        "temporarily",
        "gateway",
        "502",
        "503",
        "504",
        "net/http",
    ]
    attempt = 0
    while True:
        try:
            return run_cmd(cmd, env=env)
        except RuntimeError as exc:
            message = str(exc).lower()
            should_retry = any(token in message for token in retryable)
            if attempt >= retries or not should_retry:
                raise
            sleep_for = backoff_seconds * (2 ** attempt)
            console.print(
                f"[yellow]Falha temporaria. Nova tentativa em {sleep_for}s.[/yellow]"
            )
            time.sleep(sleep_for)
            attempt += 1


def extract_first_fna_from_zip(zip_path: Path, fasta_out_dir: Path, assembly_acc: str) -> Optional[Path]:
    import shutil
    import zipfile

    tmp_dir = zip_path.parent / f"extract_{assembly_acc}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zip_file:
        zip_file.extractall(tmp_dir)

    fna_files = list(tmp_dir.rglob("*.fna"))
    if not fna_files:
        return None

    fna = fna_files[0]
    out_path = fasta_out_dir / f"{assembly_acc}.fna"

    try:
        with fna.open("r") as handle:
            iterator = SeqIO.parse(handle, "fasta")
            _ = next(iterator, None)
    except Exception as exc:
        raise RuntimeError(f"FASTA invalido para {assembly_acc}: {fna} | erro: {exc}")

    if out_path.exists():
        out_path.unlink()

    shutil.copyfile(fna, out_path)
    return out_path


def compute_fasta_metrics(fna_path: Path) -> Dict[str, str]:
    lengths: List[int] = []
    gc_count = 0
    n_count = 0

    with fna_path.open("r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
            lengths.append(len(seq))
            gc_count += seq.count("G") + seq.count("C")
            n_count += seq.count("N")

    if not lengths:
        return {
            "observed_genome_size": "",
            "observed_contig_count": "",
            "observed_contig_n50": "",
            "observed_contig_l50": "",
            "observed_gc_percent": "",
            "observed_n_count": "",
            "observed_min_contig_len": "",
            "observed_max_contig_len": "",
        }

    total_len = sum(lengths)
    lengths_sorted = sorted(lengths, reverse=True)
    half = total_len / 2
    running = 0
    n50 = 0
    l50 = 0
    for index, length in enumerate(lengths_sorted, start=1):
        running += length
        if running >= half:
            n50 = length
            l50 = index
            break

    gc_percent = (gc_count / total_len) * 100 if total_len else 0.0

    return {
        "observed_genome_size": str(total_len),
        "observed_contig_count": str(len(lengths)),
        "observed_contig_n50": str(n50),
        "observed_contig_l50": str(l50),
        "observed_gc_percent": f"{gc_percent:.3f}",
        "observed_n_count": str(n_count),
        "observed_min_contig_len": str(min(lengths)),
        "observed_max_contig_len": str(max(lengths)),
    }


def row_from_assembly_report_line(line: Dict[str, Any], observed: Dict[str, str]) -> Row:
    assembly_info = line.get("assemblyInfo", {}) or {}
    stats = line.get("assemblyStats", {}) or {}
    organism = line.get("organism", {}) or {}

    accession = str(safe_get(line, ["accession"], "")) or str(
        safe_get(assembly_info, ["assemblyAccession"], "")
    )
    if not accession:
        accession = str(safe_get(line, ["assemblyInfo", "pairedAssembly", "accession"], ""))

    organism_name = str(safe_get(organism, ["organismName"], "")) or str(
        safe_get(organism, ["scientificName"], "")
    )
    taxid = str(safe_get(organism, ["taxId"], ""))
    bioproject = str(safe_get(assembly_info, ["bioprojectAccession"], ""))
    biosample = str(safe_get(assembly_info, ["biosample", "accession"], ""))

    submission_date = str(safe_get(assembly_info, ["submissionDate"], ""))
    release_date = str(safe_get(assembly_info, ["releaseDate"], ""))
    seq_release_date = str(safe_get(assembly_info, ["seqReleaseDate"], ""))
    seed_year = ""
    if release_date:
        seed_year = release_date[:4]
    elif submission_date:
        seed_year = submission_date[:4]
    elif seq_release_date:
        seed_year = seq_release_date[:4]

    assembly_level = str(safe_get(assembly_info, ["assemblyLevel"], ""))
    sequencing_tech = str(safe_get(assembly_info, ["sequencingTech"], ""))

    geo_loc_name = str(safe_get(assembly_info, ["biosample", "geoLocName"], ""))
    country = geo_loc_name.split(":", 1)[0].strip() if geo_loc_name else ""

    genome_size = str(safe_get(stats, ["totalSequenceLength"], ""))
    total_ungapped_length = str(safe_get(stats, ["totalUngappedLength"], ""))
    gaps_between_scaffolds = str(safe_get(stats, ["gapsBetweenScaffoldsCount"], ""))
    number_of_chromosomes = str(safe_get(stats, ["totalNumberOfChromosomes"], ""))
    number_of_organelles = str(safe_get(stats, ["numberOfOrganelles"], ""))
    number_of_scaffolds = str(safe_get(stats, ["numberOfScaffolds"], ""))
    scaffold_n50 = str(safe_get(stats, ["scaffoldN50"], ""))
    scaffold_l50 = str(safe_get(stats, ["scaffoldL50"], ""))
    number_of_contigs = str(safe_get(stats, ["numberOfContigs"], ""))
    contig_n50 = str(safe_get(stats, ["contigN50"], ""))
    contig_l50 = str(safe_get(stats, ["contigL50"], ""))
    gc_percent = str(safe_get(stats, ["gcPercent"], ""))
    source = infer_source_from_accession(accession)

    return Row(
        sequencing_tech=sequencing_tech,
        geo_loc_name=geo_loc_name,
        country=country,
        assembly_accession=accession,
        assembly_submission_date=submission_date,
        assembly_release_date=release_date,
        assembly_seq_release_date=seq_release_date,
        seed_year=seed_year,
        organism_name=organism_name,
        taxid=taxid,
        bioproject=bioproject,
        biosample=biosample,
        source=source,
        assembly_level=assembly_level,
        genome_size=genome_size,
        total_ungapped_length=total_ungapped_length,
        gaps_between_scaffolds=gaps_between_scaffolds,
        number_of_chromosomes=number_of_chromosomes,
        number_of_organelles=number_of_organelles,
        number_of_scaffolds=number_of_scaffolds,
        scaffold_n50=scaffold_n50,
        scaffold_l50=scaffold_l50,
        number_of_contigs=number_of_contigs,
        contig_n50=contig_n50,
        contig_l50=contig_l50,
        gc_percent=gc_percent,
        observed_genome_size=observed.get("observed_genome_size", ""),
        observed_contig_count=observed.get("observed_contig_count", ""),
        observed_contig_n50=observed.get("observed_contig_n50", ""),
        observed_contig_l50=observed.get("observed_contig_l50", ""),
        observed_gc_percent=observed.get("observed_gc_percent", ""),
        observed_n_count=observed.get("observed_n_count", ""),
        observed_min_contig_len=observed.get("observed_min_contig_len", ""),
        observed_max_contig_len=observed.get("observed_max_contig_len", ""),
    )


def read_assembly_report_from_extracted(zip_path: Path, assembly_acc: str) -> Dict[str, Any]:
    import zipfile

    tmp_dir = zip_path.parent / f"extract_{assembly_acc}"
    if not tmp_dir.exists():
        tmp_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path, "r") as zip_file:
            zip_file.extractall(tmp_dir)

    report_path = tmp_dir / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
    if not report_path.exists():
        raise RuntimeError(f"assembly_data_report.jsonl nao encontrado em {tmp_dir}")

    with report_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            if str(obj.get("accession", "")) == assembly_acc:
                return obj

    with report_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                return json.loads(line)

    raise RuntimeError(f"assembly_data_report.jsonl vazio para {assembly_acc}")


def write_csv(rows: List[Row], out_csv: Path) -> None:
    headers = [
        "sequencing_tech",
        "geo_loc_name",
        "country",
        "assembly_accession",
        "assembly_submission_date",
        "assembly_release_date",
        "assembly_seq_release_date",
        "seed_year",
        "organism_name",
        "taxid",
        "bioproject",
        "biosample",
        "source",
        "assembly_level",
        "genome_size",
        "total_ungapped_length",
        "gaps_between_scaffolds",
        "number_of_chromosomes",
        "number_of_organelles",
        "number_of_scaffolds",
        "scaffold_n50",
        "scaffold_l50",
        "number_of_contigs",
        "contig_n50",
        "contig_l50",
        "gc_percent",
        "observed_genome_size",
        "observed_contig_count",
        "observed_contig_n50",
        "observed_contig_l50",
        "observed_gc_percent",
        "observed_n_count",
        "observed_min_contig_len",
        "observed_max_contig_len",
    ]

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for row in rows:
            writer.writerow([
                row.sequencing_tech,
                row.geo_loc_name,
                row.country,
                row.assembly_accession,
                row.assembly_submission_date,
                row.assembly_release_date,
                row.assembly_seq_release_date,
                row.seed_year,
                row.organism_name,
                row.taxid,
                row.bioproject,
                row.biosample,
                row.source,
                row.assembly_level,
                row.genome_size,
                row.total_ungapped_length,
                row.gaps_between_scaffolds,
                row.number_of_chromosomes,
                row.number_of_organelles,
                row.number_of_scaffolds,
                row.scaffold_n50,
                row.scaffold_l50,
                row.number_of_contigs,
                row.contig_n50,
                row.contig_l50,
                row.gc_percent,
                row.observed_genome_size,
                row.observed_contig_count,
                row.observed_contig_n50,
                row.observed_contig_l50,
                row.observed_gc_percent,
                row.observed_n_count,
                row.observed_min_contig_len,
                row.observed_max_contig_len,
            ])


def write_failures(failures: List[Dict[str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["assembly_accession", "error"])
        for failure in failures:
            writer.writerow([
                failure.get("assembly_accession", ""),
                failure.get("error", ""),
            ])


def status_flag(ok: bool, label: str) -> str:
    return label if ok else "-"


def render_dashboard(state: DashboardState, overall_progress: Progress, item_progress: Progress) -> Group:
    counters = state.counters_by_source()

    summary = Table.grid(expand=True)
    summary.add_column(justify="left", ratio=1)
    summary.add_column(justify="left", ratio=1)
    summary.add_column(justify="left", ratio=1)
    summary.add_row(
        "job: oculto",
        "taxon: oculto",
        f"hits: {state.total_hits}",
    )
    summary.add_row(
        f"concluidos: {state.completed}/{state.total_hits}",
        f"falhas: {state.failures}",
        f"pendentes: {state.pending_count()}",
    )
    summary.add_row(
        f"genbank ok: {counters.get('GenBank', 0)}",
        f"refseq ok: {counters.get('RefSeq', 0)}",
        f"atual: {state.current_accession}",
    )

    status_table = Table(expand=True)
    status_table.add_column("Accession", no_wrap=True)
    status_table.add_column("Banco", no_wrap=True)
    status_table.add_column("ZIP", width=5)
    status_table.add_column("FASTA", width=7)
    status_table.add_column("Report", width=8)
    status_table.add_column("Status", no_wrap=True)
    status_table.add_column("Nota", overflow="fold")

    for item in state.statuses[-12:]:
        status_table.add_row(
            item.accession,
            item.source,
            status_flag(item.zip_ready, "OK"),
            status_flag(item.fasta_ready, "OK"),
            status_flag(item.report_ready, "OK"),
            item.state,
            item.note,
        )

    return Group(
        Panel(summary, title="1_bd :: download NCBI", border_style="cyan"),
        overall_progress,
        item_progress,
        Panel(status_table, title="Itens recentes", border_style="green"),
    )


# -----------------------------
# Main job
# -----------------------------
def main(cfg_path: str) -> None:
    cfg_file = Path(cfg_path).resolve()
    load_local_dotenv(cfg_file)
    with cfg_file.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}

    base = (cfg_file.parent / cfg["output"]["base_dir"]).resolve()
    dirs = ensure_dirs(base)
    fasta_dir = dirs["fasta"]
    meta_dir = dirs["meta"]
    pkgs_dir = dirs["pkgs"]

    _email, api_key = resolve_ncbi_credentials(cfg)
    if not command_available("datasets"):
        raise RuntimeError("Comando 'datasets' nao encontrado no ambiente atual")

    taxon_value = str(cfg["query"]["value"])
    job_name = str(cfg["request_name"])
    run_cfg = cfg.get("run", {}) or {}
    retries = int(run_cfg.get("retries", 3))
    backoff_seconds = int(run_cfg.get("retry_backoff_seconds", 5))
    strict = bool(run_cfg.get("strict", False))
    local_cfg = cfg.get("local_metrics", {}) or {}
    local_metrics = bool(local_cfg.get("enable", True))

    summary_cmd = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        taxon_value,
        "--as-json-lines",
    ]
    if api_key.strip():
        summary_cmd += ["--api-key", api_key.strip()]

    console.print("[cyan]Consultando assemblies do taxon no NCBI...[/cyan]")
    res = run_cmd_with_retries(
        summary_cmd,
        retries=retries,
        backoff_seconds=backoff_seconds,
    )

    reports = parse_jsonl_stdout(res.stdout)
    accessions: List[str] = []
    for report in reports:
        accession = safe_get(report, ["assembly", "assembly_accession"])
        if accession:
            accessions.append(str(accession))
    if not accessions:
        for report in reports:
            accession = report.get("accession")
            if accession:
                accessions.append(str(accession))
    if not accessions:
        raise RuntimeError("Nenhum assembly encontrado para o taxon informado")

    rows: List[Row] = []
    failures: List[Dict[str, str]] = []
    out_csv = meta_dir / "assemblies.csv"
    failed_csv = meta_dir / "failed_assemblies.csv"
    done_flag = meta_dir / ".job_done"

    state = DashboardState(
        job_name=job_name,
        taxon_value=taxon_value,
        total_hits=len(accessions),
        statuses=[
            ItemStatus(accession=acc, source=infer_source_from_accession(acc))
            for acc in accessions
        ],
    )

    overall_progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        TextColumn("[white]{task.fields[current]}"),
        console=console,
        expand=True,
    )
    overall_task_id = overall_progress.add_task(
        "geral",
        total=len(accessions),
        current="iniciando",
    )
    item_total = 3 if local_metrics else 2
    item_progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold green]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TextColumn("[white]{task.fields[current]}"),
        console=console,
        expand=True,
    )
    item_task_id = item_progress.add_task(
        "item",
        total=item_total,
        current="aguardando",
    )

    with Live(render_dashboard(state, overall_progress, item_progress), console=console, refresh_per_second=8) as live:
        overall_progress.start()
        item_progress.start()
        for item in state.statuses:
            accession = item.accession
            state.current_accession = accession
            overall_progress.update(overall_task_id, current=accession)
            item_progress.reset(item_task_id, total=item_total, completed=0, current="baixando pacote")
            item.state = "running"
            item.note = "baixando pacote"
            live.update(render_dashboard(state, overall_progress, item_progress))

            zip_path = pkgs_dir / f"{accession}.zip"
            download_cmd = [
                "datasets",
                "download",
                "genome",
                "accession",
                accession,
                "--include",
                "genome",
                "--filename",
                str(zip_path),
            ]
            if api_key.strip():
                download_cmd += ["--api-key", api_key.strip()]

            try:
                if not is_valid_zip(zip_path):
                    if zip_path.exists():
                        zip_path.unlink()
                    run_cmd_with_retries(
                        download_cmd,
                        retries=retries,
                        backoff_seconds=backoff_seconds,
                    )
                item.zip_ready = True
                item.note = "pacote pronto"
                item_progress.update(item_task_id, advance=1, current="pacote pronto")
                live.update(render_dashboard(state, overall_progress, item_progress))

                observed: Dict[str, str] = {}
                if local_metrics:
                    item.note = "extraindo fasta"
                    item_progress.update(item_task_id, current="extraindo fasta")
                    live.update(render_dashboard(state, overall_progress, item_progress))
                    fna_path = extract_first_fna_from_zip(zip_path, fasta_dir, accession)
                    if fna_path:
                        item.fasta_ready = True
                        observed = compute_fasta_metrics(fna_path)
                        item.note = "fasta pronto"
                    else:
                        item.note = "pacote sem fasta"
                    item_progress.update(item_task_id, advance=1, current=item.note)
                    live.update(render_dashboard(state, overall_progress, item_progress))

                item.note = "lendo metadata"
                item_progress.update(item_task_id, current="lendo metadata")
                live.update(render_dashboard(state, overall_progress, item_progress))
                report_line = read_assembly_report_from_extracted(zip_path, accession)
                item.report_ready = True
                rows.append(row_from_assembly_report_line(report_line, observed))
                item.state = "done"
                item.note = "ok"
                item_progress.update(item_task_id, advance=1, current="metadata pronta")
                state.completed += 1
            except Exception as exc:
                error = str(exc).replace("\n", " | ")
                failures.append({"assembly_accession": accession, "error": error})
                item.state = "failed"
                item.note = error[:120]
                item_progress.update(item_task_id, current="falha")
                state.failures += 1
                if strict:
                    write_csv(rows, out_csv)
                    write_failures(failures, failed_csv)
                    raise
            finally:
                write_csv(rows, out_csv)
                write_failures(failures, failed_csv)
                overall_progress.advance(overall_task_id)
                live.update(render_dashboard(state, overall_progress, item_progress))

        overall_progress.stop()
        item_progress.stop()
        if not failures:
            done_flag.parent.mkdir(parents=True, exist_ok=True)
            done_flag.write_text("OK\n", encoding="utf-8")
        live.update(render_dashboard(state, overall_progress, item_progress))

    if failures:
        console.print(f"[yellow]Falhas registradas em {failed_csv}[/yellow]")
        raise RuntimeError(f"Fetch incompleto: {len(failures)} assemblies falharam")

    console.print(f"[bold green]Concluido[/bold green] -> {out_csv}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python pipeline/run_taxon_job.py config.yaml")
        sys.exit(1)
    main(sys.argv[1])
