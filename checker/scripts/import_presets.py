"""Reproduce checked-in constants; OpenSSL is needed only for this developer script."""

import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ghs_checker import BinaryCurve, BinaryField, analyze_curve


def extract(text, label):
    lines = text.splitlines()
    i = next(i for i, line in enumerate(lines) if line.startswith(label + ':'))
    immediate = lines[i].split(':', 1)[1].strip()
    match = re.search(r'\(0x([0-9a-fA-F]+)\)', immediate)
    if match:
        return int(match.group(1), 16)
    if immediate.isdecimal():
        return int(immediate)
    value = ''
    for line in lines[i + 1:]:
        if not re.fullmatch(r'\s+[0-9a-fA-F:]+', line):
            break
        value += line.strip().replace(':', '')
    if not value:
        raise ValueError(label)
    return int(value, 16)


def main():
    out = {}
    for degree in (163, 233, 283, 409, 571):
        for family, ending in [('B', 'r2' if degree == 163 else 'r1'), ('K', 'k1')]:
            sec = f'sect{degree}{ending}'
            text = subprocess.check_output(['openssl', 'ecparam', '-name', sec, '-param_enc', 'explicit', '-text', '-noout'], text=True)
            modulus = extract(text, 'Polynomial')
            a, b, order = (extract(text, label) for label in ('A', 'B', 'Order'))
            cofactor = int(re.search(r'Cofactor:\s+(\d+)', text).group(1))
            width = (degree + 7) // 8
            data = extract(text, 'Generator (uncompressed)').to_bytes(1 + 2 * width, 'big')
            assert data[0] == 4
            gx, gy = int.from_bytes(data[1:1 + width], 'big'), int.from_bytes(data[1 + width:], 'big')
            name = f'{family}-{degree}'
            out[name] = {'name': name, 'degree': degree, 'modulus': hex(modulus), 'a': hex(a), 'b': hex(b),
                         'subgroup_order': str(order), 'cofactor': cofactor, 'curve_order': str(order * cofactor),
                         'base_point': {'x': hex(gx), 'y': hex(gy)}, 'aliases': [sec],
                         'source': 'NIST SP 800-186 / FIPS 186-4; polynomial-basis constants cross-checked with OpenSSL',
                         'source_url': 'https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186.pdf'}
            curve = BinaryCurve(BinaryField(degree, modulus), a, b)
            assert curve.contains((gx, gy)), name

    # Exact polynomial exponents from Appendix A of the March 31, 2003 preprint.
    secret = [152,143,139,136,135,133,130,125,124,122,120,119,118,117,116,114,113,112,110,109,106,105,103,102,101,99,97,96,92,91,88,87,86,85,81,78,77,76,75,73,71,69,68,67,66,63,59,58,53,51,50,49,48,46,45,44,42,38,34,33,32,31,29,27,26,24,23,22,21,20,19,18,17,16,15,14,13,12,10,7,6,4,3,2]
    public = [160,156,155,153,152,151,150,149,148,147,146,145,143,142,141,130,129,127,126,125,124,123,120,118,112,109,104,103,102,101,99,98,97,96,93,92,91,90,88,85,83,77,74,70,68,65,64,63,62,61,60,58,57,55,50,48,45,41,38,37,36,33,31,30,27,26,24,23,22,21,20,19,17,16,14,13,10,8,7,4,3,1]
    for name, exps in [('teske-secret', secret), ('teske-public', public)]:
        assert len(exps) == len(set(exps))
        b = sum(1 << e for e in exps)
        out[name] = {'name': name, 'degree': 161, 'modulus': 'x^161+x^18+1', 'a': '0', 'b': hex(b),
                     'subgroup_order': '730750818665451459101841775429946272920385056109', 'cofactor': 4,
                     'source': 'Edlyn Teske, An Elliptic Curve Trapdoor System, 2003 preprint, Appendix A',
                     'source_url': 'https://eprint.iacr.org/2003/058', 'b_exponents': exps}
        curve = BinaryCurve(BinaryField(161, (1 << 161) | (1 << 18) | 1), 0, b)
        print(name, [(d['subfield_degree'], d['m'], d['genus']) for d in analyze_curve(curve)['descents']])
    target = ROOT / 'ghs_checker' / 'data' / 'presets.json'
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(out, indent=2) + '\n')
    print('Generated', len(out), 'validated curve definitions')


if __name__ == '__main__':
    main()
