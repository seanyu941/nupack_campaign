"""ORM models for the campaign database.

The tables split into four groups:

Catalogs, written once and reused across runs:
  ``variants``   one DNA design, carrying the scan metadata that describes how
                 it differs from the parent sequence (truncation, mutations)
  ``targets``    the strand each variant is meant to bind
  ``conditions`` one point in the parameter sweep

The sweep record:
  ``runs``            one row per sweep invocation
  ``target_energies`` free energy of a target strand on its own, one row per
                      (run, target, condition) rather than per variant
  ``results``         one row per (run, variant, target, condition) cell

Selection output:
  ``selections`` and ``selection_members``

Binding free energy is a three-body quantity:

    delta_g_binding = complex_dg - cdna_dg - target_dg

The target term does not depend on the variant, so it is computed once per
(run, target, condition) and joined in. On the 25 nt scan that is 2 evaluations
instead of 109,320 copies of the same two numbers.

``delta_g_binding_kcal`` is stored on ``results`` even though it is derivable,
because it is the column every query ranks and filters on and it needs an
index. A test checks the stored value against its components.

This module is the only place the schema is defined. ``sql/schema.sql`` and
``sql/indexes.sql`` are generated from this metadata by
``scripts/dump_schema_sql.py``, and a test checks the two stay in sync.
"""

from __future__ import annotations

from datetime import datetime, timezone

from sqlalchemy import (
    CheckConstraint,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


def utcnow() -> datetime:
    """Naive UTC timestamp. SQLite has no timezone type, so tz info is dropped."""
    return datetime.now(timezone.utc).replace(tzinfo=None)


class Base(DeclarativeBase):
    pass


class Target(Base):
    """A strand the designs are meant to bind."""

    __tablename__ = "targets"

    target_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String(64), nullable=False, unique=True)
    sequence: Mapped[str] = mapped_column(Text, nullable=False)
    description: Mapped[str | None] = mapped_column(Text, nullable=True)

    def __repr__(self) -> str:
        return f"<Target {self.name} ({len(self.sequence)} nt)>"


class Variant(Base):
    """One DNA design, plus how it was derived from the parent sequence.

    The scan columns come straight from the design CSV. They are stored rather
    than parsed out of the sequence at query time so that trend queries can
    group by them and use an index.

    Nullable position and base columns mean "no mutation at this slot", which is
    how an unmutated reference row is represented.
    """

    __tablename__ = "variants"

    variant_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String(96), nullable=False, unique=True)
    sequence: Mapped[str] = mapped_column(Text, nullable=False)
    length_nt: Mapped[int] = mapped_column(Integer, nullable=False)
    gc_content: Mapped[float] = mapped_column(Float, nullable=False)

    # Truncation. trunc_total is derived on load so trend queries can group and
    # index on it without recomputing a sum on every row.
    trunc_5prime: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    trunc_3prime: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    trunc_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    # Mutation slots. A double mutation scan fills both, a single fills the
    # first, a reference row fills neither.
    position_1: Mapped[int | None] = mapped_column(Integer, nullable=True)
    original_base_1: Mapped[str | None] = mapped_column(String(1), nullable=True)
    mutated_base_1: Mapped[str | None] = mapped_column(String(1), nullable=True)
    mutation_type_1: Mapped[str] = mapped_column(String(16), nullable=False, default="none")

    position_2: Mapped[int | None] = mapped_column(Integer, nullable=True)
    original_base_2: Mapped[str | None] = mapped_column(String(1), nullable=True)
    mutated_base_2: Mapped[str | None] = mapped_column(String(1), nullable=True)
    mutation_type_2: Mapped[str] = mapped_column(String(16), nullable=False, default="none")

    num_transitions: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    num_transversions: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    # Derived on load. num_mutations counts filled slots, variant_class buckets
    # the row, and mutation_signature is a readable label such as "G0A,C1A".
    num_mutations: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    variant_class: Mapped[str] = mapped_column(String(24), nullable=False, default="wild_type")
    mutation_signature: Mapped[str] = mapped_column(String(64), nullable=False, default="")

    # True only for the fully untruncated, unmutated design. The per-truncation
    # references are found by num_mutations = 0 instead, since a scan has one of
    # those per truncation combination.
    is_reference: Mapped[bool] = mapped_column(Integer, nullable=False, default=0)

    created_at: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=utcnow)

    results: Mapped[list[Result]] = relationship(back_populates="variant")

    __table_args__ = (
        CheckConstraint("length_nt > 0", name="ck_variants_length_positive"),
        CheckConstraint("gc_content BETWEEN 0 AND 1", name="ck_variants_gc_range"),
        CheckConstraint("num_mutations BETWEEN 0 AND 2", name="ck_variants_mutation_count"),
        CheckConstraint("trunc_5prime >= 0 AND trunc_3prime >= 0", name="ck_variants_trunc"),
        # The dimensions the trend queries group by.
        Index("idx_variants_class", "variant_class"),
        Index("idx_variants_trunc", "trunc_total", "trunc_5prime", "trunc_3prime"),
        Index("idx_variants_gc", "gc_content"),
        Index("idx_variants_position_1", "position_1"),
        # Finds the unmutated reference for a given truncation in one seek,
        # which is what the ddG join needs.
        Index("idx_variants_reference_lookup", "num_mutations", "trunc_5prime", "trunc_3prime"),
    )

    def __repr__(self) -> str:
        return f"<Variant {self.name} ({self.variant_class}, {self.length_nt} nt)>"


class VariantAlignment(Base):
    """Where a design binds on a target, and how well it pairs there.

    Derived from the two sequences rather than read from the design file, so it
    stays correct if a sequence is edited. Computed once per (variant, target)
    at load time because it does not depend on the condition.

    ``n_complementary`` is the useful one for trends: it folds truncation and
    mutation into a single count of bases that actually pair, which is closer to
    what drives binding than either axis on its own.
    """

    __tablename__ = "variant_alignments"

    variant_id: Mapped[int] = mapped_column(
        ForeignKey("variants.variant_id", ondelete="CASCADE"), primary_key=True
    )
    target_id: Mapped[int] = mapped_column(
        ForeignKey("targets.target_id", ondelete="CASCADE"), primary_key=True
    )
    target_offset: Mapped[int] = mapped_column(Integer, nullable=False)
    n_complementary: Mapped[int] = mapped_column(Integer, nullable=False)
    n_mismatch: Mapped[int] = mapped_column(Integer, nullable=False)
    match_fraction: Mapped[float] = mapped_column(Float, nullable=False)

    __table_args__ = (
        CheckConstraint("n_complementary >= 0", name="ck_alignment_complementary"),
        Index("idx_alignment_paired", "target_id", "n_complementary"),
        Index("idx_alignment_offset", "target_id", "target_offset"),
    )

    def __repr__(self) -> str:
        return (
            f"<VariantAlignment variant={self.variant_id} offset={self.target_offset} "
            f"{self.n_complementary}/{self.n_complementary + self.n_mismatch} paired>"
        )


class Condition(Base):
    """One point in the sweep grid.

    The unique constraint means re-running a sweep with an overlapping grid
    reuses existing rows instead of creating near-duplicate conditions that
    would silently split the results.
    """

    __tablename__ = "conditions"

    condition_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    temperature_c: Mapped[float] = mapped_column(Float, nullable=False)
    na_molar: Mapped[float] = mapped_column(Float, nullable=False)
    mg_molar: Mapped[float] = mapped_column(Float, nullable=False)
    material: Mapped[str] = mapped_column(String(16), nullable=False)
    ensemble: Mapped[str] = mapped_column(String(16), nullable=False)

    __table_args__ = (
        UniqueConstraint(
            "temperature_c",
            "na_molar",
            "mg_molar",
            "material",
            "ensemble",
            name="uq_conditions_point",
        ),
        Index("idx_conditions_temperature", "temperature_c"),
    )

    def __repr__(self) -> str:
        return (
            f"<Condition {self.temperature_c}C Na={self.na_molar}M "
            f"Mg={self.mg_molar}M {self.material}/{self.ensemble}>"
        )


class Run(Base):
    """One invocation of the sweep, or one import of precomputed results.

    ``git_sha`` and ``grid_hash`` are recorded so a result set can be tied back
    to the code and the grid file that produced it. ``status`` is left at
    ``running`` if the process dies, which is how a stale run is spotted later.
    """

    __tablename__ = "runs"

    run_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    started_at: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=utcnow)
    finished_at: Mapped[datetime | None] = mapped_column(DateTime, nullable=True)
    git_sha: Mapped[str] = mapped_column(String(40), nullable=False)
    engine_name: Mapped[str] = mapped_column(String(32), nullable=False)
    engine_version: Mapped[str] = mapped_column(String(32), nullable=False)
    grid_hash: Mapped[str] = mapped_column(String(64), nullable=False)
    status: Mapped[str] = mapped_column(String(16), nullable=False, default="running")
    note: Mapped[str | None] = mapped_column(Text, nullable=True)

    results: Mapped[list[Result]] = relationship(back_populates="run")

    __table_args__ = (
        CheckConstraint(
            "status IN ('running', 'complete', 'failed')", name="ck_runs_status"
        ),
        Index("idx_runs_status_started", "status", "started_at"),
    )

    def __repr__(self) -> str:
        return f"<Run {self.run_id} {self.status} engine={self.engine_name}>"


class TargetEnergy(Base):
    """Free energy of a target strand on its own.

    Separated from ``results`` because it does not depend on the variant.
    Computing it here turns one evaluation per variant into one per condition.
    """

    __tablename__ = "target_energies"

    target_energy_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    run_id: Mapped[int] = mapped_column(
        ForeignKey("runs.run_id", ondelete="CASCADE"), nullable=False
    )
    target_id: Mapped[int] = mapped_column(ForeignKey("targets.target_id"), nullable=False)
    condition_id: Mapped[int] = mapped_column(
        ForeignKey("conditions.condition_id"), nullable=False
    )
    free_energy_kcal: Mapped[float] = mapped_column(Float, nullable=False)

    __table_args__ = (
        UniqueConstraint("run_id", "target_id", "condition_id", name="uq_target_energy_cell"),
    )

    def __repr__(self) -> str:
        return (
            f"<TargetEnergy run={self.run_id} target={self.target_id} "
            f"condition={self.condition_id} dG={self.free_energy_kcal:.3f}>"
        )


class Result(Base):
    """One evaluated cell of the sweep.

    The unique constraint on (run_id, variant_id, target_id, condition_id) is
    what makes a rerun a no-op: the runner asks the database which cells are
    missing and only computes those.
    """

    __tablename__ = "results"

    result_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    run_id: Mapped[int] = mapped_column(
        ForeignKey("runs.run_id", ondelete="CASCADE"), nullable=False
    )
    variant_id: Mapped[int] = mapped_column(
        ForeignKey("variants.variant_id"), nullable=False
    )
    target_id: Mapped[int] = mapped_column(ForeignKey("targets.target_id"), nullable=False)
    condition_id: Mapped[int] = mapped_column(
        ForeignKey("conditions.condition_id"), nullable=False
    )

    # The three components. cdna_dg is the design on its own, complex_dg the
    # duplex with the target, and the target term lives in target_energies.
    cdna_dg_kcal: Mapped[float] = mapped_column(Float, nullable=False)
    complex_dg_kcal: Mapped[float] = mapped_column(Float, nullable=False)
    delta_g_binding_kcal: Mapped[float] = mapped_column(Float, nullable=False)

    # Optional extras. NUPACK gives these cheaply, an imported CSV may not have
    # them, so nothing downstream is allowed to require them.
    mfe_kcal: Mapped[float | None] = mapped_column(Float, nullable=True)
    ensemble_defect: Mapped[float | None] = mapped_column(Float, nullable=True)
    compute_seconds: Mapped[float | None] = mapped_column(Float, nullable=True)

    run: Mapped[Run] = relationship(back_populates="results")
    variant: Mapped[Variant] = relationship(back_populates="results")
    condition: Mapped[Condition] = relationship()

    __table_args__ = (
        UniqueConstraint(
            "run_id", "variant_id", "target_id", "condition_id", name="uq_results_cell"
        ),
        # Ranking is ORDER BY delta_g_binding_kcal DESC within a run and
        # condition, so this index covers the sort as well as the filter.
        Index("idx_results_ranking", "run_id", "condition_id", "delta_g_binding_kcal"),
        # Covers the per-variant lookups the ddG join and the resume check make.
        Index("idx_results_run_variant", "run_id", "variant_id"),
    )

    def __repr__(self) -> str:
        return (
            f"<Result run={self.run_id} variant={self.variant_id} "
            f"condition={self.condition_id} dGb={self.delta_g_binding_kcal:.2f}>"
        )


class Selection(Base):
    """One execution of the selection query.

    ``criteria_json`` is a snapshot of the resolved parameters rather than a
    pointer to the config file, since the config file changes. ``sql_sha256``
    covers the query text for the same reason.
    """

    __tablename__ = "selections"

    selection_id: Mapped[int] = mapped_column(Integer, primary_key=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.run_id"), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=utcnow)
    criteria_json: Mapped[str] = mapped_column(Text, nullable=False)
    sql_sha256: Mapped[str] = mapped_column(String(64), nullable=False)
    n_candidates_in: Mapped[int] = mapped_column(Integer, nullable=False)
    n_selected: Mapped[int] = mapped_column(Integer, nullable=False)

    members: Mapped[list[SelectionMember]] = relationship(
        back_populates="selection", cascade="all, delete-orphan"
    )

    __table_args__ = (Index("idx_selections_run", "run_id"),)

    def __repr__(self) -> str:
        return (
            f"<Selection {self.selection_id} run={self.run_id} "
            f"{self.n_candidates_in} -> {self.n_selected}>"
        )


class SelectionMember(Base):
    """A variant that made a given shortlist, with the metrics it was picked on."""

    __tablename__ = "selection_members"

    selection_id: Mapped[int] = mapped_column(
        ForeignKey("selections.selection_id", ondelete="CASCADE"), primary_key=True
    )
    variant_id: Mapped[int] = mapped_column(
        ForeignKey("variants.variant_id"), primary_key=True
    )
    rank: Mapped[int] = mapped_column(Integer, nullable=False)
    delta_g_binding_kcal: Mapped[float] = mapped_column(Float, nullable=False)
    ddg_vs_reference_kcal: Mapped[float | None] = mapped_column(Float, nullable=True)

    selection: Mapped[Selection] = relationship(back_populates="members")
    variant: Mapped[Variant] = relationship()

    __table_args__ = (
        UniqueConstraint("selection_id", "rank", name="uq_selection_rank"),
    )

    def __repr__(self) -> str:
        return f"<SelectionMember sel={self.selection_id} rank={self.rank}>"
