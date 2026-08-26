"""Thermodynamic engines.

The sweep runner only knows about the ``ThermoEngine`` protocol, so there are
two implementations behind it:

``NupackEngine``
    The real one. Needs a NUPACK install and a licence.

``StubEngine``
    A deterministic stand-in that returns smooth functions of sequence
    composition and condition. The numbers are not physically meaningful and
    should never be used for a design decision. It exists so that the pipeline,
    the schema, the queries and the tests can run in CI and on a machine
    without a NUPACK licence.

Both return the same ``ThermoResult``, so switching between them changes the
numbers in the database and nothing else.
"""

from __future__ import annotations

import hashlib
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np


@dataclass(frozen=True, slots=True)
class ConditionPoint:
    """One sweep condition, independent of how it is stored."""

    temperature_c: float
    na_molar: float
    mg_molar: float
    material: str
    ensemble: str

    def as_key(self) -> tuple[float, float, float, str, str]:
        return (
            self.temperature_c,
            self.na_molar,
            self.mg_molar,
            self.material,
            self.ensemble,
        )


@dataclass(frozen=True, slots=True)
class ThermoResult:
    """What one binding evaluation produces.

    ``cdna_dg_kcal`` is the design on its own and ``complex_dg_kcal`` the duplex
    with the target. Binding free energy is the difference once the target's own
    energy is taken off, which the caller supplies because it does not depend on
    the variant:

        delta_g_binding = complex_dg - cdna_dg - target_dg

    ``mfe_kcal`` is the minimum free energy structure of the design alone, so
    ``cdna_dg_kcal <= mfe_kcal`` always holds. The tests assert this.
    """

    cdna_dg_kcal: float
    complex_dg_kcal: float
    delta_g_binding_kcal: float
    mfe_kcal: float | None = None
    ensemble_defect: float | None = None


@runtime_checkable
class ThermoEngine(Protocol):
    name: str
    version: str

    def evaluate_strand(self, sequence: str, condition: ConditionPoint) -> float:
        """Free energy of one strand on its own, in kcal/mol."""
        ...

    def evaluate_binding(
        self,
        sequence: str,
        target_sequence: str,
        condition: ConditionPoint,
        target_dg_kcal: float,
    ) -> ThermoResult:
        """Free energy of the duplex, and the binding term that follows."""
        ...


COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_complement(sequence: str) -> str:
    return "".join(COMPLEMENT[base] for base in reversed(sequence.upper()))


@dataclass(frozen=True, slots=True)
class Alignment:
    """Where a design sits on its target, and how well it matches there.

    ``offset`` is the 0-based start of the binding window in the target read
    5' to 3'. ``n_complementary`` counts the positions that actually pair, so
    ``n_mismatch`` is what a mutation costs in pairing terms.
    """

    offset: int
    n_complementary: int
    length: int

    @property
    def n_mismatch(self) -> int:
        return self.length - self.n_complementary

    @property
    def match_fraction(self) -> float:
        return self.n_complementary / self.length if self.length else 0.0


class TargetAligner:
    """Finds the best binding window for a design on one target.

    A design is usually much shorter than its target, so assuming the two line
    up at position 0 is wrong. This slides the design along the target and takes
    the window with the most complementary positions.

    A design pairs with its target antiparallel, so the design's reverse
    complement read 5' to 3' is what should match a target window directly. That
    turns the search into a plain comparison and lets numpy do it in one pass per
    design instead of a nested Python loop.
    """

    def __init__(self, target_sequence: str) -> None:
        self.target = target_sequence.strip().upper()
        self._encoded = np.frombuffer(self.target.encode(), dtype=np.uint8)
        self._windows: dict[int, np.ndarray] = {}
        self._cache: dict[str, Alignment] = {}

    def _window_matrix(self, length: int) -> np.ndarray:
        """Every window of the given length, as one (n_offsets, length) array."""
        if length not in self._windows:
            n_offsets = max(len(self.target) - length + 1, 1)
            self._windows[length] = np.lib.stride_tricks.sliding_window_view(
                self._encoded, min(length, len(self.target))
            )[:n_offsets]
        return self._windows[length]

    def align(self, sequence: str) -> Alignment:
        sequence = sequence.strip().upper()
        if sequence in self._cache:
            return self._cache[sequence]

        probe = np.frombuffer(reverse_complement(sequence).encode(), dtype=np.uint8)
        windows = self._window_matrix(len(sequence))

        if windows.shape[1] < len(probe):
            # Design is longer than the target, so compare on the overlap only.
            probe = probe[: windows.shape[1]]

        matches = (windows == probe).sum(axis=1)
        offset = int(matches.argmax())

        alignment = Alignment(
            offset=offset, n_complementary=int(matches[offset]), length=len(sequence)
        )
        self._cache[sequence] = alignment
        return alignment


class NoBindingSite(ValueError):
    pass


def check_binding_site(
    sequences: Sequence[str], target_sequence: str, min_fraction: float = 0.6
) -> Alignment:
    """Confirm the designs actually have somewhere to bind on this target.

    Catches the common mistake of passing the wrong target sequence, which would
    otherwise produce a full run of results that all look plausible and are all
    meaningless. Checks the best matching design rather than every one, since a
    mutation scan is meant to contain designs that bind poorly.
    """
    aligner = TargetAligner(target_sequence)
    best = max((aligner.align(s) for s in sequences), key=lambda a: a.match_fraction)

    if best.match_fraction < min_fraction:
        raise NoBindingSite(
            f"the best of {len(sequences)} designs pairs at only "
            f"{best.match_fraction:.0%} against this target. Check the "
            "--target-sequence argument is the right strand and the right region."
        )
    return best


def gc_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    if not sequence:
        return 0.0
    return sum(base in "GC" for base in sequence) / len(sequence)


class StubEngine:
    """Deterministic placeholder for NUPACK.

    Values come from a hash of (sequence, condition), so the same input always
    gives the same output and a rerun of the sweep is byte-for-byte comparable.
    The trends are loosely sensible (more GC means a more negative dG, higher
    temperature and lower salt push it up) which is enough for the selection
    query to have something with structure to sort, but the magnitudes are made
    up.
    """

    name = "stub"
    version = "1.0"

    def __init__(self, seed: str = "nupack-campaign") -> None:
        self.seed = seed
        self._aligners: dict[str, TargetAligner] = {}

    def _aligner(self, target_sequence: str) -> TargetAligner:
        """One aligner per target, so the window matrix is built once."""
        if target_sequence not in self._aligners:
            self._aligners[target_sequence] = TargetAligner(target_sequence)
        return self._aligners[target_sequence]

    def _uniform(self, *parts: object) -> float:
        """Reproducible float in [0, 1) from the given parts."""
        payload = "|".join([self.seed, *(str(p) for p in parts)]).encode()
        digest = hashlib.blake2b(payload, digest_size=8).digest()
        return int.from_bytes(digest, "big") / 2**64

    def _strand_dg(self, sequence: str, condition: ConditionPoint) -> float:
        """Smooth, deterministic stand-in for the free energy of one strand."""
        gc = gc_fraction(sequence)
        stacks = max(len(sequence) - 1, 1)

        # Rough per-stack contribution, GC weighted more heavily than AT.
        per_stack = 0.72 + 0.90 * gc
        dg_ref = -per_stack * stacks * 0.42

        # Structure gets less stable as temperature rises above the 37 C
        # reference.
        dg_temperature = 0.045 * (condition.temperature_c - 37.0) * math.sqrt(stacks)

        # Divalent cations stabilise more per mole than monovalent ones. The
        # factor of 120 is the usual rule of thumb for a sodium equivalent.
        sodium_equivalent = condition.na_molar + 120.0 * condition.mg_molar
        dg_salt = 0.90 * math.log10(1.0 / max(sodium_equivalent, 1e-3))

        jitter = (self._uniform(sequence, condition.as_key()) - 0.5) * 1.2
        return dg_ref + dg_temperature + dg_salt + jitter

    def evaluate_strand(self, sequence: str, condition: ConditionPoint) -> float:
        return round(self._strand_dg(sequence, condition), 4)

    def evaluate_binding(
        self,
        sequence: str,
        target_sequence: str,
        condition: ConditionPoint,
        target_dg_kcal: float,
    ) -> ThermoResult:
        cdna_dg = self._strand_dg(sequence, condition)

        # The duplex is modelled as both strands plus a hybridisation term that
        # scales with how many bases actually pair at the design's best window on
        # the target. Truncating the design or mismatching a base costs pairing,
        # which is the trend the scan is looking for.
        alignment = self._aligner(target_sequence).align(sequence)
        overlap = alignment.length
        match_fraction = alignment.match_fraction

        hybridisation = -1.35 * overlap * (0.25 + 0.75 * match_fraction)
        hybridisation += 0.030 * (condition.temperature_c - 37.0) * overlap * 0.1
        hybridisation += (self._uniform("hyb", sequence, condition.as_key()) - 0.5) * 0.4

        complex_dg = cdna_dg + target_dg_kcal + hybridisation
        delta_g_binding = complex_dg - cdna_dg - target_dg_kcal

        mfe = cdna_dg + 0.05 + 0.85 * self._uniform("mfe", sequence, condition.as_key())

        defect = 0.06 + 0.30 * abs(gc_fraction(sequence) - 0.5)
        defect += 0.004 * max(condition.temperature_c - 37.0, 0.0)
        defect += 0.05 * self._uniform("defect", sequence, condition.as_key())

        return ThermoResult(
            cdna_dg_kcal=round(cdna_dg, 4),
            complex_dg_kcal=round(complex_dg, 4),
            delta_g_binding_kcal=round(delta_g_binding, 4),
            mfe_kcal=round(mfe, 4),
            ensemble_defect=round(min(max(defect, 0.005), 0.95), 4),
        )


class NupackEngine:
    """Adapter over NUPACK 4.

    Written against the 4.0 API. The import is deferred to ``__init__`` so that
    importing this module does not require NUPACK to be installed.
    """

    name = "nupack"

    def __init__(self) -> None:
        try:
            import nupack
        except ImportError as exc:  # pragma: no cover - depends on local install
            raise ImportError(
                "NUPACK is not installed. It is not on PyPI and needs a licence, "
                "see https://nupack.org. Use the stub engine to run the pipeline "
                "without it."
            ) from exc

        self._nupack = nupack
        self.version = getattr(nupack, "__version__", "unknown")
        self._model_cache: dict[tuple, object] = {}

    def _model(self, condition: ConditionPoint):  # pragma: no cover - needs NUPACK
        key = condition.as_key()
        if key not in self._model_cache:
            self._model_cache[key] = self._nupack.Model(
                material=condition.material,
                celsius=condition.temperature_c,
                sodium=condition.na_molar,
                magnesium=condition.mg_molar,
                ensemble=condition.ensemble,
            )
        return self._model_cache[key]

    def evaluate_strand(  # pragma: no cover - needs NUPACK
        self, sequence: str, condition: ConditionPoint
    ) -> float:
        # pfunc returns (partition function, free energy of the ensemble).
        _, free_energy = self._nupack.pfunc(strands=[sequence], model=self._model(condition))
        return float(free_energy)

    def evaluate_binding(  # pragma: no cover - needs NUPACK
        self,
        sequence: str,
        target_sequence: str,
        condition: ConditionPoint,
        target_dg_kcal: float,
    ) -> ThermoResult:
        nupack = self._nupack
        model = self._model(condition)

        _, cdna_dg = nupack.pfunc(strands=[sequence], model=model)
        _, complex_dg = nupack.pfunc(strands=[sequence, target_sequence], model=model)

        mfe_structures = nupack.mfe(strands=[sequence], model=model)
        mfe_energy = float(mfe_structures[0].energy)
        structure = str(mfe_structures[0].structure)
        raw_defect = float(nupack.defect(strands=[sequence], structure=structure, model=model))

        cdna_dg = float(cdna_dg)
        complex_dg = float(complex_dg)

        return ThermoResult(
            cdna_dg_kcal=cdna_dg,
            complex_dg_kcal=complex_dg,
            delta_g_binding_kcal=complex_dg - cdna_dg - target_dg_kcal,
            mfe_kcal=mfe_energy,
            # Normalised by length so the value is comparable across designs.
            ensemble_defect=raw_defect / max(len(sequence), 1),
        )


ENGINES: dict[str, type] = {"stub": StubEngine, "nupack": NupackEngine}


def get_engine(name: str) -> ThermoEngine:
    try:
        factory = ENGINES[name]
    except KeyError:
        options = ", ".join(sorted(ENGINES))
        raise ValueError(f"Unknown engine {name!r}, pick one of: {options}") from None
    return factory()
