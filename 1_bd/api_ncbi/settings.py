import os
import shutil
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import yaml
from dotenv import load_dotenv


PROJECT_MARKERS = {'.env.example', 'environment.yml', 'Snakefile'}
PLACEHOLDER_VALUES = {
    '',
    'you@example.org',
    'replace-with-your-ncbi-api-key',
}


def _find_project_root(start: Path) -> Path:
    resolved = start.resolve()
    if resolved.is_file():
        resolved = resolved.parent

    for candidate in [resolved, *resolved.parents]:
        if candidate.name == '1_bd':
            return candidate
        if any((candidate / marker).exists() for marker in PROJECT_MARKERS):
            return candidate
    return resolved


def project_root_from(start: Path) -> Path:
    return _find_project_root(start)


def read_yaml_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open('r', encoding='utf-8') as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f'Config invalido em {path}')
    return data


def write_yaml_config(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)


def ensure_local_env(start: Path) -> Tuple[Optional[Path], bool]:
    project_root = _find_project_root(start)
    env_path = project_root / '.env'
    example_path = project_root / '.env.example'

    if env_path.exists():
        return env_path, False

    if example_path.exists():
        shutil.copyfile(example_path, env_path)
        return env_path, True

    return None, False


def load_local_dotenv(start: Path) -> Optional[Path]:
    project_root = _find_project_root(start)
    env_path = project_root / '.env'
    if env_path.exists():
        load_dotenv(env_path, override=False)
        return env_path
    return None


def update_env_values(env_path: Path, values: Dict[str, str]) -> None:
    existing_lines = []
    if env_path.exists():
        existing_lines = env_path.read_text(encoding='utf-8').splitlines()

    updated_keys = set()
    output_lines = []
    for line in existing_lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#') or '=' not in line:
            output_lines.append(line)
            continue

        key, _value = line.split('=', 1)
        key = key.strip()
        if key in values:
            output_lines.append(f"{key}={values[key]}")
            updated_keys.add(key)
        else:
            output_lines.append(line)

    for key, value in values.items():
        if key not in updated_keys:
            output_lines.append(f"{key}={value}")

    env_path.write_text('\n'.join(output_lines).rstrip() + '\n', encoding='utf-8')


def command_available(command: str) -> bool:
    return shutil.which(command) is not None


def nested_get(data: Dict[str, Any], *path: str, default: str = '') -> str:
    current: Any = data
    for key in path:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
    if current is None:
        return default
    return str(current)


def sanitize_secret_value(value: str) -> str:
    normalized = value.strip()
    if normalized in PLACEHOLDER_VALUES:
        return ''
    return normalized


def resolve_ncbi_credentials(cfg: Dict[str, Any]) -> Tuple[str, str]:
    email = sanitize_secret_value(os.getenv('NCBI_EMAIL', '')) or sanitize_secret_value(
        nested_get(cfg, 'ncbi', 'email', default='')
    )
    api_key = sanitize_secret_value(os.getenv('NCBI_API_KEY', '')) or sanitize_secret_value(
        nested_get(cfg, 'ncbi', 'api_key', default='')
    )
    return email, api_key


def mask_secret(value: str) -> str:
    if not value:
        return 'nao definido'
    if len(value) <= 8:
        return '*' * len(value)
    return f'{value[:4]}...{value[-4:]}'
