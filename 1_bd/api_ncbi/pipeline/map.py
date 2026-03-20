import argparse
from pathlib import Path

import pandas as pd
import plotly.express as px


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gera mapa de bolhas por pais/ano a partir de assemblies_tratados.csv.",
    )
    parser.add_argument(
        "--input",
        default="outputs/metadata/assemblies_tratados.csv",
        help="CSV de entrada (padrao: outputs/metadata/assemblies_tratados.csv).",
    )
    parser.add_argument(
        "--output",
        default="outputs/metadata/mapa_projetos.html",
        help="HTML de saida (padrao: outputs/metadata/mapa_projetos.html).",
    )
    return parser.parse_args()


def resolve_path(base_dir: Path, path_str: str) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path
    return base_dir / path


def build_project_map(input_csv: Path, output_html: Path) -> None:
    df = pd.read_csv(input_csv)

    if "country" not in df.columns:
        raise SystemExit("Coluna 'country' nao encontrada no CSV.")

    if "geo_loc_name" in df.columns:
        mask = df["country"].isna() | (df["country"].astype(str).str.strip() == "")
        geo = df.loc[mask, "geo_loc_name"].fillna("").map(str)
        df.loc[mask, "country"] = geo.map(lambda value: value.split(":", 1)[0].strip())

    year = pd.to_numeric(df.get("seed_year"), errors="coerce")
    if "assembly_release_date" in df.columns:
        rel = pd.to_datetime(df["assembly_release_date"], errors="coerce")
        year = year.fillna(rel.dt.year)
    df["year"] = year

    df = df[df["country"].notna() & df["year"].notna()]

    yearly = (
        df.groupby(["country", "year"], dropna=True)
        .size()
        .reset_index(name="project_count")
    )
    years = sorted(yearly["year"].dropna().astype(int).unique().tolist())
    if not years:
        raise SystemExit("Nenhum ano valido encontrado para gerar o mapa.")

    all_countries = sorted(yearly["country"].dropna().astype(str).unique().tolist())
    full_index = pd.MultiIndex.from_product([all_countries, years], names=["country", "year"])
    agg = yearly.set_index(["country", "year"]).reindex(full_index).reset_index()
    agg["project_count"] = agg["project_count"].fillna(0).astype(int)
    agg["project_count"] = agg.groupby("country")["project_count"].cumsum()

    fig = px.scatter_geo(
        agg,
        locations="country",
        locationmode="country names",
        size="project_count",
        color="project_count",
        hover_name="country",
        hover_data={"project_count": True, "year": True},
        projection="natural earth",
        size_max=40,
        title="Projetos acumulados por pais",
        animation_frame="year",
    )

    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html))
    print(f"OK: {output_html}")
    print("Abra o HTML no navegador para explorar o mapa animado.")


def main() -> None:
    args = parse_args()
    base_dir = Path(__file__).resolve().parents[1]
    input_csv = resolve_path(base_dir, args.input)
    output_html = resolve_path(base_dir, args.output)
    build_project_map(input_csv, output_html)


if __name__ == "__main__":
    main()
