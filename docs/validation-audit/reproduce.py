import hashlib,json,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from importlib.metadata import version
from benchmarks.checkpoints import batch_path,read_checkpoint
from benchmarks.expanded import summarize,generate,evaluate,_seed
from benchmarks.report import calibration

out=Path('benchmark_results/validation')
manifest=json.loads((out/'config.json').read_text()); c=manifest['config']
report=dict(audited_at=datetime.now(timezone.utc).isoformat(),configuration=c,errors=[],checks={},replays=[])
report['checks']['source_hashes_match']={p:Path(p).exists() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest for p,digest in manifest['source_hashes'].items()}
report['checks']['dependency_versions_match']={p:version(p)==v for p,v in manifest['versions'].items()}
report['checks']['historical_runner_fingerprint_present']=(out/'runner_metadata.json').exists()
records=[];terms=[];file_hashes={};n=0
for scenario in c['scenarios']:
    shard=out/'shards'/scenario
    expected_manifest=dict(manifest,config=dict(c,scenarios=[scenario]))
    assert json.loads((shard/'config.json').read_text())==expected_manifest,scenario+' manifest mismatch'
    expected={batch_path(shard,scenario,e,start) for e in c['odds_ratios'] for start in range(0,c['replicates'],c['batch_size'])}
    assert expected==set(shard.glob('batch_*.json')),scenario+' checkpoint names mismatch'
    for e in c['odds_ratios']:
        for start in range(0,c['replicates'],c['batch_size']):
            p=batch_path(shard,scenario,e,start)
            payload=read_checkpoint(p,c,scenario,e,start)
            digest=hashlib.sha256(p.read_bytes()).hexdigest()
            assert digest==hashlib.sha256((out/p.name).read_bytes()).hexdigest(),p.name+' root differs'
            file_hashes[p.name]=digest
            records.extend(payload['replicates']);terms.extend(payload['terms']);n+=1
    print('Validated',scenario,n,'batches',flush=True)
assert set(file_hashes)=={p.name for p in out.glob('batch_*.json')}
report['checks']['validated_batches']=n
raw=pd.DataFrame(records).sort_values(['scenario','odds_ratio','replicate','bins','method']).reset_index(drop=True)
table=pd.DataFrame(terms).sort_values(['scenario','odds_ratio','replicate','bins','Term_ID']).reset_index(drop=True)
assert not raw.duplicated(['scenario','odds_ratio','replicate','bins','method']).any()
assert not table.duplicated(['scenario','odds_ratio','replicate','bins','Term_ID']).any()
report['checks']['metric_rows']=len(raw);report['checks']['term_rows']=len(table)
for name,frame,keys in [('replicates',raw,['scenario','odds_ratio','replicate','bins','method']),('terms',table,['scenario','odds_ratio','replicate','bins','Term_ID']),('summary',summarize(raw,c['seed']),['scenario','odds_ratio','bins','method']),('calibration',calibration(table),['scenario','odds_ratio','bins','method','threshold'])]:
    saved=pd.read_csv(out/(name+'.tsv'),sep='\t',low_memory=False)
    for col in ['conditional_null','marginal_null']:
        if col in frame:
            frame[col]=frame[col].astype('boolean')
            saved[col]=saved[col].astype('boolean')
    try:
        pd.testing.assert_frame_equal(frame.sort_values(keys).reset_index(drop=True),saved.sort_values(keys).reset_index(drop=True),check_dtype=False,check_exact=False,rtol=1e-8,atol=1e-10)
        report['checks'][name+'_matches_reconstruction']=True
    except AssertionError as error:
        report['checks'][name+'_matches_reconstruction']=False
        report['errors'].append(dict(table=name,error=str(error)[:1500]))
        destination=Path('docs/validation-audit');destination.mkdir(exist_ok=True)
        frame.to_csv(destination/(name+'_reconstructed.tsv'),sep='\t',index=False)
    print('Reconstructed',name,flush=True)
for scenario in c['scenarios']:
    for e in c['odds_ratios']:
        data=generate(np.random.default_rng(_seed(c['seed'],scenario,e,0,0)),scenario,c['n_families'],e)
        actual,fraction=evaluate(data,2,c['permutations'],_seed(c['seed'],scenario,e,0,1,2))
        saved=table[(table.scenario==scenario)&(table.odds_ratio==e)&(table.replicate==0)&(table.bins==2)].set_index('Term_ID')
        for column in actual.columns:
            left=actual[column].sort_index();right=saved[column].sort_index()
            if column in ['conditional_null','marginal_null']:
                left=left.astype('boolean');right=right.astype('boolean')
            pd.testing.assert_series_equal(left,right,check_dtype=False,check_names=False,check_exact=False,rtol=1e-8,atol=1e-10)
        report['replays'].append(dict(scenario=scenario,odds_ratio=e,replicate=0,bins=2,matches=True))
    print('Replayed',scenario,flush=True)
report['checks']['sampled_replays_match']=len(report['replays'])
report['historical_provenance']='Consult per-file source checks: two package files differ from the recorded versions. Dependency versions and the numerical modules used by the benchmark match. Historical orchestration fingerprint is absent; structural consistency and sampled replay do not prove provenance of all historical tasks.'
destination=Path('docs/validation-audit');destination.mkdir(exist_ok=True)
(destination/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
(destination/'checkpoint_hashes.json').write_text(json.dumps(file_hashes,indent=2,sort_keys=True)+'\n')
pd.read_csv(out/'summary.tsv',sep='\t').to_csv(destination/'summary.tsv',sep='\t',index=False)
print('AUDIT COMPLETED',report['checks'],flush=True)
