"""Runs inside the CIBERSORTx container (Python 3.6). Copied from Evin's
deconvolution_LM22_20261006/run_02_bmode_recovery/workflow/native.py (11.03); the only change is that the HiRes
label is an argument. Credentials are read from a mounted JSON file and replaced by [redacted] in every log line."""
import argparse
import datetime
import json
import subprocess
import sys
import threading

p = argparse.ArgumentParser()
p.add_argument('--stage', choices=['fractions', 'hires'], required=True)
p.add_argument('--auth', required=True)
p.add_argument('--threads', type=int, required=True)
p.add_argument('--label', default='GTEx_LM22')
p.add_argument('--window', type=int, default=88)   # HiRes window (samples); 88 as in the GTEx run, smaller for GSE60424
a = p.parse_args()
with open(a.auth) as stream:
    auth = json.load(stream)
secrets = [auth['username'], auth['token']]
if not all(isinstance(x, str) and x.strip() for x in secrets):
    raise ValueError('Account credentials are empty.')


def emit(line):
    for value in sorted(secrets, key=len, reverse=True):
        line = line.replace(value, '[redacted]')
    print(line, flush=True)


# same flags as the earlier GTEx run: B-mode batch correction on, quantile normalisation off
command = ['/src/CIBERSORTxFractions' if a.stage == 'fractions' else '/src/CIBERSORTxHiRes',
           '--username', auth['username'], '--token', auth['token'], '--mixture', 'mixture.tsv',
           '--sigmatrix', 'signature.tsv', '--QN', 'FALSE', '--rmbatchBmode', 'TRUE',
           '--outdir', '/src/outdir']
if a.stage == 'fractions':
    command += ['--perm', '100', '--absolute', 'FALSE']
else:
    command += ['--cibresults', 'fractions.tsv', '--subsetgenes', 'gene_subset.txt',
                '--threads', str(a.threads), '--nsampling', '100', '--nsampling2', '10',
                '--window', str(a.window), '--variableonly', 'FALSE', '--label', a.label]
stopped = threading.Event()


def heartbeat():
    while not stopped.wait(60):
        emit('[%s] Native %s remains active.' % (datetime.datetime.utcnow().isoformat() + 'Z', a.stage))


threading.Thread(target=heartbeat, daemon=True).start()
emit('Start native ' + a.stage + '; credentials hidden; B-mode TRUE; QN FALSE.')
try:
    process = subprocess.Popen(command, cwd='/src', stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True, bufsize=1)
    for line in process.stdout:
        emit(line.rstrip('\r\n'))
    status = process.wait()
finally:
    stopped.set()
emit('Native exit status: ' + str(status))
sys.exit(status)
