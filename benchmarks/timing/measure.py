import json, time, platform, hashlib
from pathlib import Path
from importlib.metadata import version
import numpy as np
from benchmarks.expanded import generate, evaluate, _seed
config = json.loads(Path('benchmarks/configs/validation.json').read_text())
rows=[]
started=time.perf_counter()
for scenario in config['scenarios']:
    data=generate(np.random.default_rng(_seed(config['seed'], scenario, 1, 0, 0)), scenario, config['n_families'], 1)
    for bins in config['bins']:
        tick=time.perf_counter()
        table, fraction=evaluate(data,bins,config['permutations'],_seed(config['seed'],scenario,1,0,1,bins))
        row=dict(scenario=scenario,bins=bins,seconds=time.perf_counter()-tick,informative_strata_fraction=fraction)
        rows.append(row)
        print(json.dumps(row),flush=True)
        Path('benchmarks/timing/measurements.json').write_text(json.dumps(rows,indent=2)+'\n')
result=dict(config=config,measurement_effect=1,replicate=0,measurements=rows,elapsed_seconds=time.perf_counter()-started,estimated_compute_hours=sum(r['seconds'] for r in rows)*config['replicates']*len(config['odds_ratios'])/3600,python=platform.python_version(),platform=platform.platform(),versions={p:version(p) for p in ['numpy','pandas','scipy']},source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path('benchmarks/expanded.py'),*sorted(Path('scbenrich').rglob('*.py'))]})
Path('benchmarks/timing/runtime.json').write_text(json.dumps(result,indent=2)+'\n')
print('ESTIMATED_HOURS',result['estimated_compute_hours'],flush=True)
