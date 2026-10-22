"""Run metadata that remains readable after failed analyses."""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
from importlib.metadata import version


def build_manifest(args, inputs=None, outputs=None, status='running'):
    root = Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root, capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=root, capture_output=True, text=True, check=True).stdout)
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    return {'status': status, 'arguments': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
            'python_version': platform.python_version(), 'software_version': version('scbenrich'),
            'dependencies': {name: version(name) for name in ('numpy', 'pandas', 'scipy')},
            'git_commit': commit, 'git_dirty': dirty, 'rng': 'numpy.default_rng/PCG64',
            'code_hashes': {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted(root.rglob('*.py'))},
            'inputs': inputs or {}, 'outputs': outputs or []}


def write_manifest(path, manifest):
    path = Path(path)
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)
