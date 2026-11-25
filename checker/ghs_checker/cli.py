"""Dependency-free command-line interface."""

import argparse
import json
from pathlib import Path
import sys

from .analysis import analyze_curve, field_bounds
from .curve import BinaryCurve
from .field import BinaryField, integer, parse_element
from .isogeny import scan_isogenies
from .presets import get_preset, list_presets


def load_curve(data):
    if not isinstance(data, dict):
        raise ValueError("Curve input must be a JSON object")
    if data.get('characteristic', 2) != 2:
        raise ValueError("This checker supports only characteristic two")
    if data.get('model', 'binary-ordinary') not in ('binary-ordinary', 'y^2+x*y=x^3+a*x^2+b'):
        raise ValueError("Only the ordinary binary model y^2+x*y=x^3+a*x^2+b is supported")
    for key in ('degree', 'modulus', 'a', 'b'):
        if key not in data:
            raise ValueError(f"Missing curve input: {key}")
    field = BinaryField(integer(data['degree'], 'Field degree'), parse_element(data['modulus']))
    curve = BinaryCurve(field, field.element(data['a']), field.element(data['b']))
    order = parse_element(data['subgroup_order']) if data.get('subgroup_order') is not None else None
    if order is not None and order < 2:
        raise ValueError("Subgroup order must be >=2")
    if data.get('curve_order') is not None:
        full_order = parse_element(data['curve_order'])
        if order is not None and full_order % order:
            raise ValueError("Claimed subgroup order does not divide the supplied curve order")
    if data.get('base_point') is not None:
        point = data['base_point']
        if not isinstance(point, dict) or not {'x', 'y'} <= point.keys():
            raise ValueError("base_point must contain x and y")
        point = field.element(point['x']), field.element(point['y'])
        if not curve.contains(point):
            raise ValueError("Supplied base point does not lie on the curve")
    metadata = {key: data[key] for key in ('name', 'source', 'source_url', 'cofactor', 'curve_order') if key in data}
    if data.get('base_point') is not None:
        metadata['base_point_on_curve_verified'] = True
        metadata['base_point_order_verified'] = False
    return curve, order, metadata


def input_curve(args):
    if args.preset:
        if any(getattr(args, key, None) is not None for key in ('degree', 'modulus', 'a', 'b', 'subgroup_order')):
            raise ValueError("Use a preset unchanged, or supply custom parameters without --preset")
        data = get_preset(args.preset)
    elif args.input:
        if any(getattr(args, key, None) is not None for key in ('degree', 'modulus', 'a', 'b', 'subgroup_order')):
            raise ValueError("Use --input unchanged, or supply custom parameters without --input")
        path = Path(args.input)
        if path.stat().st_size > 1024 * 1024:
            raise ValueError("Curve input must be at most 1 MiB")
        data = json.loads(path.read_text())
    else:
        data = {key: getattr(args, key) for key in ('degree', 'modulus', 'a', 'b') if getattr(args, key) is not None}
        data['subgroup_order'] = args.subgroup_order
    return load_curve(data)


def format_report(report):
    if 'field' not in report:
        return json.dumps(report, indent=2)
    name = report['metadata'].get('name', 'custom curve')
    lines = [f"{name} over GF(2^{report['field']['degree']})", 'Classical GHS screening; no security verdict.', '']
    nodes = report.get('nodes', [{'index': 0, 'depth': 0, 'descents': report['descents']}])
    for node in nodes:
        if 'nodes' in report:
            lines.append(f"Curve {node['index']} | path length {node['depth']}")
        for d in node['descents']:
            genus = str(d['genus']) if d['genus'] is not None else 'not asserted'
            if len(genus) > 32:
                genus = f"{d['genus_formula']} ({len(str(d['genus']))} decimal digits)"
            lines.append(f"  {d['base_field']}: n={d['extension_degree']}, m={d['m']}, genus={genus}")
            lines.append(f"    {d['screening']}; subgroup: {d['subgroup_capacity']['status']}")
    if not report['descents']:
        lines.append('  No proper subfield exists.')
    if 'search' in report:
        search = report['search']
        lines.extend(['', f"Checked {search['curves_checked']} curves using degree-3 isogenies only.",
                      f"Candidate curve indices: {report['candidate_nodes']}",
                      f"Curve budget reached: {search['curve_budget_reached']}; depth limit reached: {search['depth_limit_reached']}"])
    lines.extend(['', 'Descent map: not constructed. Target subgroup image: not verified.',
                  'DLP cost: not estimated. Arbitrary covers and generalized GHS: not assessed.'])
    return '\n'.join(lines)


def parser():
    root = argparse.ArgumentParser(description='Classical GHS screening and bounded binary 3-isogeny exploration')
    root.add_argument('--version', action='version', version='ghs-checker 0.1.0')
    commands = root.add_subparsers(dest='command', required=True)
    listing = commands.add_parser('list-presets', help='List NIST and published-example presets')
    listing.add_argument('--json', action='store_true')
    bounds = commands.add_parser('field', help='Subfields and prime-degree classical GHS bounds')
    bounds.add_argument('--degree', type=int, required=True)
    bounds.add_argument('--json', action='store_true')
    for command in ('check', 'scan'):
        cmd = commands.add_parser(command, help='Screen a curve' if command == 'check' else 'Screen a bounded degree-3 isogeny neighborhood')
        source = cmd.add_mutually_exclusive_group()
        source.add_argument('--preset')
        source.add_argument('--input', help='Path to a curve JSON file')
        cmd.add_argument('--degree', type=int)
        cmd.add_argument('--modulus', help='Polynomial or integer encoding')
        cmd.add_argument('--a')
        cmd.add_argument('--b')
        cmd.add_argument('--subgroup-order')
        cmd.add_argument('--subfield-degree', type=int, action='append', help='Repeatable; defaults to all proper subfields')
        cmd.add_argument('--max-genus', type=int, default=32, help='Screening threshold only; default 32')
        cmd.add_argument('--json', action='store_true')
        if command == 'scan':
            cmd.add_argument('--max-curves', type=int, default=8)
            cmd.add_argument('--max-depth', type=int, default=3)
    return root


def main(argv=None):
    root = parser()
    args = root.parse_args(argv)
    try:
        if args.command == 'list-presets':
            presets = list_presets()
            report = {name: {key: item[key] for key in ('degree', 'source', 'source_url')} for name, item in presets.items()}
            print(json.dumps(report, indent=2) if args.json else '\n'.join(f"{name:16} GF(2^{item['degree']})" for name, item in presets.items()))
            return 0
        if args.command == 'field':
            report = field_bounds(args.degree)
        else:
            curve, order, metadata = input_curve(args)
            if args.max_genus < 1:
                raise ValueError('--max-genus must be positive')
            if args.command == 'check':
                report = analyze_curve(curve, order, args.subfield_degree, args.max_genus, metadata)
            else:
                report = scan_isogenies(curve, order, args.subfield_degree, args.max_genus,
                                        args.max_curves, args.max_depth, metadata)
        print(json.dumps(report, indent=2) if args.json else format_report(report))
        return 0
    except (ValueError, OSError, ZeroDivisionError, ArithmeticError) as exc:
        message = f'{type(exc).__name__}: {exc}'
        if args.json:
            print(json.dumps({'error': message, 'security_certification': False}), file=sys.stderr)
        else:
            print(f'ghs-checker: {message}', file=sys.stderr)
        return 2
