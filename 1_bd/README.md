<p align="center">
  <img src="https://img.shields.io/badge/versao-12.1-blue" alt="Versao">
  <img src="https://img.shields.io/badge/status-em%20desenvolvimento-yellow" alt="Status">
  <img src="https://img.shields.io/badge/feito%20com-Python%203.11-blue?logo=python&logoColor=white" alt="Python 3.11">
  <img src="https://img.shields.io/badge/workflow-Snakemake-brightgreen" alt="Snakemake">
  <img src="https://img.shields.io/badge/ambiente-Conda-44A833?logo=anaconda&logoColor=white" alt="Conda">
  <img src="https://img.shields.io/badge/dados-NCBI%20Datasets-0B6E4F" alt="NCBI Datasets">
  <img src="https://img.shields.io/badge/licenca-MIT-green" alt="Licenca">
</p>

# BIOINFO101: Genome download wrappers for NCBI

Pipeline do `seq_pipe` para busca de assemblies no NCBI, extracao de metadados, deduplicacao e geracao de visualizacoes para o modulo `1_bd`.

## Stack de desenvolvimento

- Python 3.11
- Conda/Miniconda para ambiente reproduzivel
- Snakemake para orquestracao do workflow
- NCBI Datasets CLI para consulta e download de assemblies
- Rich para interface de terminal e acompanhamento de progresso
- Pandas para tratamento tabular
- Biopython para leitura e metricas de FASTA
- Matplotlib para graficos
- Plotly para visualizacao interativa em HTML
- YAML + dotenv para configuracao segura e separacao de segredos

## Sumario

- [O que o fluxo faz](#o-que-o-fluxo-faz)
- [Arquitetura do modulo](#arquitetura-do-modulo)
- [Estrutura principal](#estrutura-principal)
- [Credenciais e seguranca](#credenciais-e-seguranca)
- [Ambiente Conda](#ambiente-conda)
- [Como rodar](#como-rodar)
- [Entradas do fluxo](#entradas-do-fluxo)
- [Saidas do fluxo](#saidas-do-fluxo)
- [Como o pipeline se unifica](#como-o-pipeline-se-unifica)
- [Observacoes operacionais](#observacoes-operacionais)

## O que o fluxo faz

1. Consulta o NCBI para um taxon.
2. Baixa os assemblies encontrados com `datasets`.
3. Extrai FASTA e metadados locais.
4. Gera `assemblies.csv` com a tabela bruta.
5. Deduplica os registros e gera `assemblies_tratados.csv`.
6. Gera graficos de sobreposicao pre-deduplicacao.
7. Gera mapa HTML de projetos por pais/ano.

## Arquitetura do modulo

A organizacao do `1_bd` segue esta estrutura:

```text
1_bd/
|-- README.md
|-- Snakefile
|-- environment.yml
|-- .env.example
`-- api_ncbi/
    |-- api_ncbi.py
    |-- config.yaml
    |-- settings.py
    |-- processing.py
    |-- pipeline/
    |   |-- run_taxon_job.py
    |   |-- analysis.py
    |   |-- dedup_graph.py
    |   `-- map.py
    `-- outputs/
        |-- metadata/
        |-- fasta/
        `-- packages/
```

Papel de cada parte:

- `Snakefile`: define o DAG do workflow no Snakemake.
- `environment.yml`: define o ambiente Conda reproduzivel.
- `.env.example`: documenta as variaveis locais esperadas.
- `api_ncbi/api_ncbi.py`: launcher interativo que encadeia o fluxo completo.
- `api_ncbi/config.yaml`: parametros versionaveis do job.
- `api_ncbi/settings.py`: leitura de `.env`, config e resolucao de paths.
- `api_ncbi/processing.py`: regras compartilhadas de deduplicacao.
- `api_ncbi/pipeline/`: scripts das etapas operacionais.
- `api_ncbi/outputs/`: resultados e artefatos locais gerados pela execucao.

## Estrutura principal

- `Snakefile`: orquestracao do fluxo no Snakemake.
- `environment.yml`: ambiente Conda do `1_bd`.
- `.env.example`: modelo de variaveis sensiveis locais.
- `README.md`: documentacao do modulo.
- `api_ncbi/api_ncbi.py`: launcher interativo do fluxo completo.
- `api_ncbi/config.yaml`: configuracao segura e versionavel do job.
- `api_ncbi/settings.py`: leitura de config, `.env` e verificacoes de ambiente.
- `api_ncbi/processing.py`: logica compartilhada de deduplicacao.
- `api_ncbi/pipeline/run_taxon_job.py`: etapa de download e montagem da tabela bruta.
- `api_ncbi/pipeline/analysis.py`: etapa de deduplicacao.
- `api_ncbi/pipeline/dedup_graph.py`: etapa de graficos.
- `api_ncbi/pipeline/map.py`: etapa de mapa.

## Credenciais e seguranca

As credenciais nao devem ir para o Git.

Arquivo esperado localmente:

```env
NCBI_EMAIL=seu_email@example.org
NCBI_API_KEY=sua_chave_ncbi
```

O launcher automatiza a parte inicial do `.env`:

1. se `1_bd/.env` nao existir, ele cria automaticamente a partir de `.env.example`
2. pede email e API key na execucao
3. pode salvar as credenciais localmente em `1_bd/.env`
4. mantem os segredos fora do `config.yaml`

O `.gitignore` na raiz de `seq_pipe` ja ignora `.env`, `outputs/`, `.snakemake/` e `__pycache__/`.

## Ambiente Conda

Criar o ambiente:

```bash
conda env create -f environment.yml
```

Ativar o ambiente:

```bash
conda activate seq-pipe-1bd
```

Pacotes principais:

- `ncbi-datasets-cli`
- `snakemake-minimal`
- `rich`
- `biopython`
- `pandas`
- `plotly`
- `matplotlib`
- `python-dotenv`

## Como rodar

### Opcao 1: launcher interativo

Roda o fluxo completo a partir de um unico arquivo:

```bash
python api_ncbi/api_ncbi.py
```

Esse launcher:

1. garante que exista um `.env` local
2. le `.env`
3. pede ou reutiliza credenciais e parametros salvos
4. grava `config.yaml` sem segredos
5. chama as 4 etapas do `1_bd` em sequencia

### Opcao 2: Snakemake

Roda o mesmo fluxo de forma declarativa:

```bash
snakemake -s Snakefile --cores 1
```

Dry-run:

```bash
snakemake -n -s Snakefile
```

## Entradas do fluxo

Entradas principais:

- `api_ncbi/config.yaml`
- `.env`
- taxon informado pelo usuario
- comando externo `datasets`

## Saidas do fluxo

Saidas principais em `api_ncbi/outputs/metadata/`:

- `assemblies.csv`
- `failed_assemblies.csv`
- `assemblies_tratados.csv`
- `dedup_summary.csv`
- `dedup_story_before.png`
- `dedup_upset_before.png`
- `mapa_projetos.html`
- `.job_done`

Saidas auxiliares locais:

- `api_ncbi/outputs/fasta/*.fna`
- `api_ncbi/outputs/packages/*.zip`

## Como o pipeline se unifica

Arquivo que unifica tudo em modo interativo:

- `api_ncbi/api_ncbi.py`

Ele encadeia estas etapas:

1. `pipeline/run_taxon_job.py`
2. `pipeline/analysis.py`
3. `pipeline/dedup_graph.py`
4. `pipeline/map.py`

No modo Snakemake, a unificacao acontece pelo:

- `Snakefile`

Ou seja:

- `api_ncbi.py` unifica o fluxo em execucao interativa e manual.
- `Snakefile` unifica o mesmo fluxo em execucao reproduzivel por regras.

<p align="center">
  <img src="https://visitor-badge.laobi.icu/badge?page_id=seq_pipe.1_bd" alt="Visitantes">
</p>