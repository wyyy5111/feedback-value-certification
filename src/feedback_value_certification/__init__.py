"""Rational value certification for contractive feedback models."""

from .certificates import certificate, projected_certificate
from .scalar import EvaluationLimitError, Model, evaluate_exact

__version__ = "1.0.0"
__all__ = [
    "certificate",
    "projected_certificate",
    "Model",
    "evaluate_exact",
    "EvaluationLimitError",
]
