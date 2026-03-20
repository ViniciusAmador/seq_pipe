import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

from rich.console import Console
from rich.panel import Panel
from rich.prompt import Confirm, Prompt
from rich.table import Table

from settings import (
    command_available,
    ensure_local_env,
    load_local_dotenv,
    mask_secret,
    project_root_from,
    read_yaml_config,
    resolve_ncbi_credentials,
    update_env_values,
    write_yaml_config,
)

console = Console()


def run_python_script(
    script: Path, args: List[str], cwd: Path, env: Optional[Dict[str, str]] = None
) -> None:
    subprocess.run(
        [sys.executable, str(script), *args],
        check=True,
        cwd=str(cwd),
        env=env,
    )


def build_safe_config(
    existing_cfg: Dict[str, object],
    taxon: str,
    request_name: str,
    base_dir: str,
    enable_local: bool,
) -> Dict[str, object]:
    run_cfg = existing_cfg.get("run", {}) if isinstance(existing_cfg, dict) else {}
    if not isinstance(run_cfg, dict):
        run_cfg = {}

    return {
        "query": {
            "mode": "taxon",
            "value": str(taxon),
        },
        "request_name": request_name,
        "output": {
            "base_dir": base_dir,
        },
        "run": {
            "retries": int(run_cfg.get("retries", 3)),
            "retry_backoff_seconds": int(run_cfg.get("retry_backoff_seconds", 5)),
            "strict": bool(run_cfg.get("strict", False)),
        },
        "local_metrics": {
            "enable": enable_local,
        },
    }


def ask_saved_value(label: str, current: str, *, secret: bool = False) -> str:
    if current:
        if Confirm.ask(f"Usar valor salvo para {label.lower()}?", default=True):
            return current
    return Prompt.ask(label, default="", password=secret)


def private_status(value: str) -> str:
    return "definido" if value.strip() else "nao definido"


def main() -> None:
    console.rule("[bold cyan]1_bd :: API NCBI[/bold cyan]")

    config_path = Path(__file__).with_name("config.yaml")
    project_root = project_root_from(config_path)
    env_path, env_created = ensure_local_env(config_path)
    load_local_dotenv(config_path)
    existing_cfg = read_yaml_config(config_path)

    if env_created and env_path is not None:
        console.print(
            f"[yellow].env criado automaticamente a partir de {project_root / '.env.example'}[/yellow]"
        )

    env_email, env_api_key = resolve_ncbi_credentials(existing_cfg)
    email = ask_saved_value("Email do NCBI (opcional)", env_email)
    api_key = env_api_key
    if env_api_key:
        if not Confirm.ask("Usar API key salva no .env?", default=True):
            api_key = Prompt.ask(
                "Informe sua API key do NCBI (opcional)",
                default="",
                password=True,
            )
    else:
        api_key = Prompt.ask(
            "Informe sua API key do NCBI (opcional)",
            default="",
            password=True,
        )

    if env_path is not None and (email.strip() or api_key.strip()):
        save_default = env_created or not env_email or not env_api_key
        if Confirm.ask("Salvar credenciais localmente em 1_bd/.env?", default=save_default):
            update_env_values(
                env_path,
                {
                    "NCBI_EMAIL": email.strip(),
                    "NCBI_API_KEY": api_key.strip(),
                },
            )
            console.print(f"[green]Credenciais salvas em {env_path}[/green]")

    query_cfg = existing_cfg.get("query", {}) if isinstance(existing_cfg, dict) else {}
    output_cfg = existing_cfg.get("output", {}) if isinstance(existing_cfg, dict) else {}
    local_cfg = existing_cfg.get("local_metrics", {}) if isinstance(existing_cfg, dict) else {}

    taxon_default = str(query_cfg.get("value", "Paracidovorax"))
    request_default = str(existing_cfg.get("request_name", "paracidovorax"))
    output_default = str(output_cfg.get("base_dir", "outputs"))
    local_default = bool(local_cfg.get("enable", True))

    taxon = ask_saved_value("Taxon (TaxID ou nome cientifico)", taxon_default)
    request_name = ask_saved_value("Nome do job", request_default)
    base_dir = Prompt.ask("Diretorio de saida", default=output_default)
    enable_local = Confirm.ask(
        "Calcular metricas locais a partir do FASTA?",
        default=local_default,
    )

    if not command_available("datasets"):
        raise SystemExit("Comando 'datasets' nao encontrado. Ative o ambiente Conda de 1_bd antes de rodar.")

    cfg = build_safe_config(existing_cfg, taxon, request_name, base_dir, enable_local)
    write_yaml_config(config_path, cfg)

    summary = Table.grid(padding=(0, 2))
    summary.add_row("Job", "oculto")
    summary.add_row("Taxon", "oculto")
    summary.add_row("Saida", base_dir)
    summary.add_row("Email", private_status(email))
    summary.add_row("API key", private_status(api_key))
    summary.add_row(".env local", str(env_path) if env_path else "nao disponivel")
    summary.add_row("Config segura", str(config_path))
    summary.add_row("Segredos no arquivo", "nao")
    console.print(Panel(summary, title="Resumo da execucao", border_style="cyan"))

    runtime_env = os.environ.copy()
    runtime_env["NCBI_EMAIL"] = email
    if api_key.strip():
        runtime_env["NCBI_API_KEY"] = api_key.strip()

    work_dir = config_path.parent
    fetch_script = Path(__file__).with_name("pipeline").joinpath("run_taxon_job.py")
    analysis_script = Path(__file__).with_name("pipeline").joinpath("analysis.py")
    graph_script = Path(__file__).with_name("pipeline").joinpath("dedup_graph.py")
    map_script = Path(__file__).with_name("pipeline").joinpath("map.py")
    metadata_dir = Path(base_dir) / "metadata"
    raw_csv = metadata_dir / "assemblies.csv"
    treated_csv = metadata_dir / "assemblies_tratados.csv"
    map_html = metadata_dir / "mapa_projetos.html"

    console.print("[cyan]Etapa 1/4: download, extracao e tabela bruta[/cyan]")
    run_python_script(fetch_script, [str(config_path)], work_dir, runtime_env)

    console.print("[cyan]Etapa 2/4: deduplicacao[/cyan]")
    run_python_script(
        analysis_script,
        ["--input", str(raw_csv), "--output", str(treated_csv)],
        work_dir,
        runtime_env,
    )

    console.print("[cyan]Etapa 3/4: graficos pre-deduplicacao[/cyan]")
    run_python_script(
        graph_script,
        ["--input", str(raw_csv), "--outdir", str(metadata_dir)],
        work_dir,
        runtime_env,
    )

    console.print("[cyan]Etapa 4/4: mapa de projetos[/cyan]")
    run_python_script(
        map_script,
        ["--input", str(treated_csv), "--output", str(map_html)],
        work_dir,
        runtime_env,
    )

    console.print("[bold green]Fluxo 1_bd concluido[/bold green]")


if __name__ == "__main__":
    main()
