#!/usr/bin/env python3
"""Regenerate all catalogues needed for SR(9)<=109; check integer weights.

Uses Python's standard library, a C++17 compiler and OpenMP. No optimization
solver or floating-point computation is used in the verification.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent

def read(path):
    raw=path.read_bytes()
    if len(raw)%16:raise ValueError('truncated catalogue')
    rows=[int.from_bytes(raw[i:i+16],'little') for i in range(0,len(raw),16)]
    if len(rows)!=len(set(rows)):raise ValueError('duplicate subsets')
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--threads',type=int,default=24)
    p.add_argument('--work',type=Path,default=HERE/'verification')
    p.add_argument('--lower-only',action='store_true')
    a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True)
    lower={}
    lower_file=HERE/'lower_partitions.json'
    if lower_file.exists():
        for key,w in json.loads(lower_file.read_text()).items():
            n,k=w['n'],w['k'];parts=w['parts']
            if len(parts)>k or any(not part for part in parts):raise ValueError('invalid number of parts')
            if sorted(x for part in parts for x in part)!=list(range(1,n+1)):raise ValueError('invalid partition')
            for part in parts:
                sums=[x+y for x,y in itertools.combinations_with_replacement(part,2)]
                diffs=[y-x for x,y in itertools.combinations(sorted(part),2)]
                if len(sums)!=len(set(sums)) or len(diffs)!=len(set(diffs)):raise ValueError('non-Sidon part')
            canonical=sorted((sorted(part) for part in parts),key=lambda part:part[0])
            digest=hashlib.sha256(json.dumps(canonical,separators=(',',':')).encode()).hexdigest()
            if digest!=w['partition_sha256']:raise ValueError('partition hash mismatch')
            lower[key]=dict(n=n,k=k,lower_bound=n+1,partition_sha256=w['partition_sha256'])
            print(f'VERIFIED SR({k}) >= {n+1}',flush=True)
        (a.work/'lower-report.json').write_text(json.dumps(lower,indent=2)+'\n')
    if a.lower_only:return
    start=time.monotonic()
    data=json.loads((HERE/'certificates.json').read_text())
    binary=a.work/'enumerate_sums'
    subprocess.run(['g++','-O3','-std=c++17','-fopenmp',str(HERE/'enumerate_sums.cpp'),'-o',str(binary)],check=True)
    def generate(size,blocked,name):
        path=a.work/(name+'.bin')
        if blocked is None:hi,lo='full','0'
        else:hi,lo=f'{blocked>>64:x}',f'{blocked&((1<<64)-1):x}'
        subprocess.run([str(binary.resolve()),'109',str(size),str(a.threads),hi,lo,str(path)],check=True)
        return read(path)
    large=generate(13,None,'size13');larger=generate(14,None,'size14')
    expected={label:sum(1<<(x-1) for x in xs) for label,xs in data['large_sets'].items()}
    for i in range(3):
        if expected[f'A_{i}']!=expected['A_0']<<i or expected[f'B_{i}']!=expected['B_0']<<i:
            raise ValueError('translation mismatch')
    if set(large)!=set(expected.values()) or larger:raise ValueError('incorrect large-set catalogue')
    pairs={frozenset((x,y)) for x,y in itertools.combinations(expected,2) if not expected[x]&expected[y]}
    if pairs!={frozenset(('A_0','B_1')),frozenset(('A_1','B_2'))}:raise ValueError('incorrect disjoint pairs')
    results=[]
    for c in data['certificates']:
        blocked=0
        for mask in c['blocked_masks']:blocked|=int(mask)
        rows=generate(12,blocked,'case'+str(c['source_case']))
        if len(rows)!=c['count']:raise ValueError('incorrect subset count')
        digest=hashlib.sha256(b''.join(r.to_bytes(16,'little') for r in sorted(rows))).hexdigest()
        if digest!=c['sorted_sha256']:raise ValueError('catalogue mismatch')
        w=c['weights'];den=c['denominator']
        if len(w)!=109 or den<=0 or any(type(x)!=int or x<0 for x in w):raise ValueError('invalid weights')
        minimum=None
        for row in rows:
            if row&blocked or row>>109:raise ValueError('invalid subset')
            # Recheck pair sums explicitly, independently of the generator's
            # bit operations, before checking its weight inequality.
            xs=[x+1 for x in range(109) if row>>x&1]
            sums=[x+y for x,y in itertools.combinations_with_replacement(xs,2)]
            if len(xs)!=12 or len(sums)!=len(set(sums)):raise ValueError('non-Sidon row')
            weight=sum(w[x-1] for x in xs)
            if weight<den:raise ValueError('row weight too small')
            minimum=weight if minimum is None else min(minimum,weight)
        if sum(w)!=c['total_weight'] or sum(w)>=c['goal']*den:raise ValueError('insufficient total weight')
        result=dict(label=c['label'],rows=len(rows),minimum=minimum,total=sum(w),goal=c['goal'])
        results.append(result);print('VERIFIED',json.dumps(result),flush=True)
    # Reflection x -> 110-x maps the other four cases to these four.
    reflect=lambda mask:sum(1<<(108-x) for x in range(109) if mask>>x&1)
    for i in range(3):
        if reflect(expected[f'A_{i}'])!=expected[f'B_{2-i}']:raise ValueError('reflection mismatch')
    report=dict(status='verified',bound='SR(9) <= 109',seconds=time.monotonic()-start,
                threads=a.threads,size13_count=len(large),size14_count=0,certificates=results,lower=lower)
    (a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('VERIFIED SR(9) <= 109',flush=True)

if __name__=='__main__':main()
