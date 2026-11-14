"""One CIBERSORTx run (Fractions or HiRes) on any mixture, with the pinned containers Evin used.
Based on Evin's deconvolution_LM22_20261006/run_02_bmode_recovery/workflow/compute.py (11.03), without its
GTEx-only checks (670 samples, 22 cell types), so it also runs per batch group and on GSE60424.

  fractions: --mixture <genes x samples TSV, first column GeneSymbol> --signature LM22.tsv --outdir <new dir>
             -> <outdir>/fractions.tsv (B-mode adjusted fractions, one row per sample)
  hires:     same + --fractions <fractions TSV> --subset <gene list>  -> <outdir>/native/*_Window88.txt

The mixture, signature and gene subset are copied into <outdir>/inputs/ (signature restricted to the genes in the
mixture, cell-type names made file-name safe, as compute.py does). The credentials file must be readable only by
its owner; it is mounted read-only and never written to a log."""
import argparse
import csv
import json
import re
import shutil
import subprocess
from pathlib import Path

CONTAINERS = Path('/oak/stanford/groups/smontgom/epadhi/AoU/cell_type_deconvolution_benchmark/containers')
SIF = {'fractions': CONTAINERS / 'cibersortx.fractions.9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc.sif',
       'hires': CONTAINERS / 'cibersortx.hires.e8da6850311d163e33a343d29a0d2ffc8b18c1ec4604995b8285c7bc2017c83e.sif'}
SINGULARITY = '/usr/bin/singularity'
HELPER = Path(__file__).with_name('cibersortx_native.py').resolve()


def safe_name(x):
    # same conversion as compute.py: anything other than letters, digits and _ becomes _
    return re.sub(r'_+', '_', re.sub(r'[^A-Za-z0-9_]', '_', x)).strip('_')


def prepare(mixture, signature, subset, fractions, inputs):
    inputs.mkdir(parents=True)
    with mixture.open() as stream:
        header = stream.readline().rstrip('\n').split('\t')
        genes = {line.split('\t', 1)[0] for line in stream}
    if header[0] != 'GeneSymbol':
        raise ValueError('The mixture must start with a GeneSymbol column.')
    shutil.copyfile(mixture, inputs / 'mixture.tsv')
    with signature.open() as stream:
        rows = list(csv.reader(stream, delimiter='\t'))
    cells = [safe_name(x) for x in rows[0][1:]]
    shared = [row for row in rows[1:] if row[0] in genes]
    with (inputs / 'signature.tsv').open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(['GeneSymbol'] + cells)
        writer.writerows(shared)
    if subset is not None:
        shutil.copyfile(subset, inputs / 'gene_subset.txt')
    report = {'mixture': str(mixture), 'samples': len(header) - 1, 'mixture_genes': len(genes),
              'signature': str(signature), 'cell_types': cells, 'signature_genes_in_mixture': len(shared)}
    (inputs / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print('prepared: %d samples, %d genes, %d LM22 genes in the mixture' % (report['samples'], len(genes), len(shared)),
          flush=True)
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['fractions', 'hires'], required=True)
    p.add_argument('--mixture', type=Path, required=True)
    p.add_argument('--signature', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--auth', type=Path, required=True)
    p.add_argument('--subset', type=Path)
    p.add_argument('--fractions', type=Path)
    p.add_argument('--label', default='GTEx_LM22')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--window', type=int, default=88)
    a = p.parse_args()
    if a.outdir.exists():
        raise ValueError('Output directory already exists; use a new one: ' + str(a.outdir))
    if not a.auth.is_file() or a.auth.stat().st_mode & 0o077:
        raise ValueError('Credentials file is missing or readable by others.')
    if a.stage == 'hires' and (a.subset is None or a.fractions is None):
        raise ValueError('HiRes needs --subset and --fractions.')

    report = prepare(a.mixture.resolve(), a.signature.resolve(), a.subset, a.fractions, a.outdir / 'inputs')
    application, native = a.outdir / 'application', a.outdir / 'native'
    application.mkdir()
    native.mkdir()
    if a.stage == 'hires':
        shutil.copyfile(a.fractions, native / 'fractions.tsv')

    # the container's /src is read-only: copy the application out first, then run it from a writable copy
    subprocess.run([SINGULARITY, 'exec', '--cleanenv', '--containall', '--bind',
                    str(application.resolve()) + ':/run/application', str(SIF[a.stage]), '/bin/sh', '-c',
                    'cp -R /src/. /run/application/'], check=True)
    command = [SINGULARITY, 'exec', '--cleanenv', '--containall', '--pwd', '/src',
               '--bind', str(application.resolve()) + ':/src',
               '--bind', str((a.outdir / 'inputs').resolve()) + ':/src/data:ro',
               '--bind', str(native.resolve()) + ':/src/outdir',
               '--bind', str(a.auth.resolve()) + ':/run/auth.json:ro',
               '--bind', str(HELPER) + ':/run/native.py:ro',
               '--env', 'OMP_NUM_THREADS=%d,OPENBLAS_NUM_THREADS=%d,MKL_NUM_THREADS=%d' % ((a.threads,) * 3),
               str(SIF[a.stage]), 'python3', '/run/native.py', '--stage', a.stage,
               '--auth', '/run/auth.json', '--threads', str(a.threads), '--label', a.label,
               '--window', str(a.window)]
    print('start CIBERSORTx %s (B-mode on, QN off)' % a.stage, flush=True)
    with (a.outdir / 'native.log').open('w') as stream:
        subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
    log = (a.outdir / 'native.log').read_text()
    if 'Native exit status: 0' not in log or any(x in log for x in ['Abort!', 'Execution halted', 'Assertion failed']):
        raise ValueError('CIBERSORTx reported an error; see ' + str(a.outdir / 'native.log'))

    if a.stage == 'fractions':
        # with B-mode on, the fraction result is CIBERSORTx_Adjusted.txt (Mixtures_Adjusted holds genes, not fractions)
        result = native / 'CIBERSORTx_Adjusted.txt'
        if not result.is_file():
            raise ValueError('CIBERSORTx_Adjusted.txt is missing.')
        shutil.copyfile(result, a.outdir / 'fractions.tsv')
        print('fractions written: ' + str((a.outdir / 'fractions.tsv').resolve()), flush=True)
    else:
        n = len(list(native.glob('*_Window*.txt')))
        if n != len(report['cell_types']):
            raise ValueError('Expected %d HiRes matrices, found %d.' % (len(report['cell_types']), n))
        print('HiRes matrices written: %d in %s' % (n, native.resolve()), flush=True)


if __name__ == '__main__':
    main()
