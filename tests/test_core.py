from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess
import sys
import unittest

from feedback_value_certification import (
    EvaluationLimitError,
    Model,
    evaluate_exact,
    projected_certificate,
)

ROOT = Path(__file__).resolve().parents[1]


def example(name):
    return json.loads((ROOT / "examples" / name).read_text())


class CertificateTests(unittest.TestCase):
    def test_rps_grouped_certificate(self):
        raw = example("rps_certificate.json")
        answer = projected_certificate(raw["kernels"], raw["projections"], raw["contraction"])
        self.assertEqual(answer["status"], "CERTIFIED_UNIFORM_BLOCK_WORD_DECAY")
        self.assertEqual(answer["block_length"], 1)
        self.assertEqual(answer["block_survival_upper"], "179/180")
        self.assertTrue(
            all(
                sorted(map(len, g["collision_groups"])) == [1, 2]
                for g in answer["projection_geometry"]
            )
        )

    def test_unnormalized_row_is_rejected(self):
        raw = example("rps_certificate.json")
        raw["kernels"]["0"]["0"][0][0] = "2"
        with self.assertRaises(ValueError):
            projected_certificate(raw["kernels"], raw["projections"], raw["contraction"])

    def test_float_probability_is_rejected(self):
        raw = example("rps_certificate.json")
        raw["kernels"]["0"]["0"][0][0] = 0.5
        with self.assertRaises(TypeError):
            projected_certificate(raw["kernels"], raw["projections"], raw["contraction"])

    def test_empty_kernel_is_rejected(self):
        with self.assertRaises(ValueError):
            projected_certificate({}, [{"0": "0"}], "2/3")


class ScalarTests(unittest.TestCase):
    def setUp(self):
        self.raw = example("reactive_binary.json")

    def test_saved_short_horizon_value(self):
        answer = evaluate_exact(Model(self.raw, "tracking"), 6)
        self.assertEqual(answer["value"], "7467729/9765625")
        self.assertEqual(answer["interval"], [answer["value"]] * 2)

    def test_actual_tag_at_threshold(self):
        model = Model(self.raw, "tracking")
        self.assertEqual(model.tag(Fraction(1, 2)), 1)
        self.assertEqual(evaluate_exact(model, 0)["value"], "0")
        self.assertEqual(evaluate_exact(model, 1)["value"], "0")

    def test_policy_has_its_own_feedback(self):
        answer = evaluate_exact(Model(self.raw, "contrarian"), 6)
        self.assertEqual(answer["value"], "-2048649/9765625")

    def test_construction_limit_is_reported(self):
        with self.assertRaises(EvaluationLimitError):
            evaluate_exact(Model(self.raw, "tracking"), 6, max_atoms=1)

    def test_invalid_distribution_is_rejected(self):
        self.raw["initial_state_law"] = ["2", "-1"]
        with self.assertRaises(ValueError):
            Model(self.raw, "tracking")

    def test_horizon_must_be_a_nonnegative_integer(self):
        model = Model(self.raw, "tracking")
        for value in [-1, 0.5, True]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                evaluate_exact(model, value)


class CommandLineTests(unittest.TestCase):
    def test_scalar_module_entry_point(self):
        call = subprocess.run(
            [
                sys.executable,
                "-m",
                "feedback_value_certification",
                "scalar",
                "--input",
                str(ROOT / "examples/reactive_binary.json"),
                "--policy",
                "tracking",
                "--horizon",
                "6",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(json.loads(call.stdout)["value"], "7467729/9765625")

    def test_module_help(self):
        call = subprocess.run(
            [sys.executable, "-m", "feedback_value_certification", "--help"],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("certificate", call.stdout)


@unittest.skipUnless(
    os.environ.get("RPS_VALUE_EXECUTABLE"), "Build the C++ executable to check RPS"
)
class RpsTests(unittest.TestCase):
    def run_rps(self, *args):
        return subprocess.run(
            [os.environ["RPS_VALUE_EXECUTABLE"], *args], capture_output=True, text=True, timeout=20
        )

    def test_saved_short_horizon_enclosures(self):
        call = self.run_rps("--T", "6", "--G", "8")
        self.assertEqual(call.returncode, 0, call.stderr)
        answer = json.loads(call.stdout)
        for name, bounds in {
            "mode": (-1024781209657, -1024781209650),
            "BR": (-889402099912, -889402099905),
        }.items():
            self.assertEqual(
                (answer[name]["lower_numerator"], answer[name]["upper_numerator"]), bounds
            )
        self.assertEqual(answer["exact_prefix_rounds"], 6)
        self.assertEqual(answer["grid_layers"], 0)

    def test_causal_continuation_and_difference(self):
        call = self.run_rps("--T", "9", "--G", "8")
        self.assertEqual(call.returncode, 0, call.stderr)
        answer = json.loads(call.stdout)
        self.assertEqual(answer["exact_prefix_rounds"], 8)
        self.assertEqual(answer["grid_layers"], 1)
        self.assertTrue(answer["same_target_piece_for_all_next_WL_states"])
        q = answer["fixed_point_Q"]
        for key in ["mode", "BR"]:
            lower, upper = answer[key]["lower_numerator"], answer[key]["upper_numerator"]
            self.assertLessEqual(-9 * q, lower)
            self.assertLessEqual(lower, upper)
            self.assertLessEqual(upper, 9 * q)
        self.assertEqual(
            answer["BR_minus_mode"]["lower_numerator"],
            answer["BR"]["lower_numerator"] - answer["mode"]["upper_numerator"],
        )

    def test_memory_preflight_stops_before_allocation(self):
        call = self.run_rps("--T", "6", "--G", "100000")
        self.assertEqual(call.returncode, 2)
        self.assertIn("RESOURCE", json.loads(call.stdout)["status"])

    def test_parameters_are_explicit(self):
        call = self.run_rps("--T", "6")
        self.assertNotEqual(call.returncode, 0)


if __name__ == "__main__":
    unittest.main()
