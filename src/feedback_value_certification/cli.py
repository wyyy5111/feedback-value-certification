"""Command-line access to certificates and exact scalar values."""

import argparse
import json
from pathlib import Path

from . import __version__
from .certificates import projected_certificate
from .scalar import EvaluationLimitError, Model, evaluate_exact


def _positive_integer(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("Expected a positive integer")
    return number


def _horizon(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("Expected a nonnegative integer")
    return number


def _reject_float(value):
    raise ValueError('Write non-integer rational inputs as strings, for example "2/3"')


def load_input(path):
    with Path(path).open(encoding="utf-8") as stream:
        data = json.load(stream, parse_float=_reject_float, parse_constant=_reject_float)
    if not isinstance(data, dict):
        raise ValueError("The input must be a JSON object")
    return data


def main(argv=None):
    parser = argparse.ArgumentParser(description="Exact rational feedback-model evaluation")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    cert = commands.add_parser("certificate", help="Check projected separation and block survival")
    cert.add_argument("--input", "-i", required=True, type=Path)
    scalar = commands.add_parser("scalar", help="Compute an exact scalar expected reward")
    scalar.add_argument("--input", "-i", required=True, type=Path)
    scalar.add_argument("--policy", required=True)
    scalar.add_argument("--horizon", "-T", required=True, type=_horizon)
    scalar.add_argument("--max-cuts", type=_positive_integer, default=12000)
    scalar.add_argument("--max-atoms", type=_positive_integer, default=200000)
    args = parser.parse_args(argv)
    try:
        raw = load_input(args.input)
        if args.command == "certificate":
            result = projected_certificate(raw["kernels"], raw["projections"], raw["contraction"])
        else:
            result = evaluate_exact(
                Model(raw, args.policy),
                args.horizon,
                max_cuts=args.max_cuts,
                max_atoms=args.max_atoms,
            )
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except EvaluationLimitError as error:
        print(json.dumps({"status": "RESOURCE_LIMIT", "message": str(error)}))
        return 2
    except (ValueError, TypeError, KeyError, OSError, ZeroDivisionError, IndexError) as error:
        print(json.dumps({"status": "INVALID_INPUT", "message": str(error)}))
        return 1
