import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from processing import apply_dedup_pipeline, load_csv_rows, write_csv_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Deduplica assemblies por biosample e remove linhas sem sequencing_tech.",
    )
    parser.add_argument(
        "--input",
        default="outputs/metadata/assemblies.csv",
        help="Caminho do CSV de entrada (padrao: outputs/metadata/assemblies.csv).",
    )
    parser.add_argument(
        "--output",
        default="outputs/metadata/assemblies_tratados.csv",
        help="Caminho do CSV de saida (padrao: outputs/metadata/assemblies_tratados.csv).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base_dir = Path(__file__).resolve().parents[1]
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.is_absolute():
        input_path = base_dir / input_path
    if not output_path.is_absolute():
        output_path = base_dir / output_path

    rows = load_csv_rows(input_path)
    if not rows:
        raise SystemExit(f"Nenhuma linha encontrada em {input_path}")

    fieldnames = list(rows[0].keys())
    treated, _stats = apply_dedup_pipeline(rows)
    write_csv_rows(output_path, treated, fieldnames)
    print(f"OK: {output_path} (linhas: {len(treated)})")


if __name__ == "__main__":
    main()
