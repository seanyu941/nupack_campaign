"""Engine and session helpers.

Everything else in the package takes a ``Session`` and does not care where the
database lives, so the same code runs against the committed SQLite file and
against the in-memory database the tests use.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from sqlalchemy import Engine, create_engine, event, text
from sqlalchemy.dialects import sqlite as sqlite_dialect
from sqlalchemy.orm import Session, sessionmaker
from sqlalchemy.schema import CreateIndex, CreateTable

from .models import Base

DEFAULT_DB_URL = "sqlite:///data/campaign.db"

# Window functions in sql/select_candidates.sql need SQLite 3.25 or newer.
MIN_SQLITE_VERSION = (3, 25, 0)


class UnsupportedSQLiteError(RuntimeError):
    pass


def check_sqlite_version() -> None:
    version = tuple(int(part) for part in sqlite3.sqlite_version.split("."))
    if version < MIN_SQLITE_VERSION:
        wanted = ".".join(str(p) for p in MIN_SQLITE_VERSION)
        raise UnsupportedSQLiteError(
            f"SQLite {sqlite3.sqlite_version} is too old, the selection query uses "
            f"window functions and needs {wanted} or newer."
        )


def make_engine(url: str = DEFAULT_DB_URL, echo: bool = False) -> Engine:
    """Build an Engine and turn on the SQLite settings this project assumes.

    Foreign keys are off by default in SQLite, so the ON DELETE CASCADE on
    results and selection_members would be ignored without the pragma.
    """
    if url.startswith("sqlite:///") and ":memory:" not in url:
        Path(url.removeprefix("sqlite:///")).parent.mkdir(parents=True, exist_ok=True)

    engine = create_engine(url, echo=echo, future=True)

    if engine.dialect.name == "sqlite":
        check_sqlite_version()

        @event.listens_for(engine, "connect")
        def _set_sqlite_pragmas(dbapi_connection, _record):  # noqa: ANN001
            cursor = dbapi_connection.cursor()
            cursor.execute("PRAGMA foreign_keys=ON")
            # WAL keeps the notebook readable while a sweep is still writing.
            cursor.execute("PRAGMA journal_mode=WAL")
            cursor.execute("PRAGMA synchronous=NORMAL")
            cursor.close()

    return engine


def init_db(engine: Engine) -> None:
    """Create tables and indexes if they are not there yet."""
    Base.metadata.create_all(engine)


def make_session_factory(engine: Engine) -> sessionmaker[Session]:
    return sessionmaker(bind=engine, future=True, expire_on_commit=False)


@contextmanager
def session_scope(engine: Engine) -> Iterator[Session]:
    """Commit on the way out, roll back if the body raises."""
    session = make_session_factory(engine)()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def _tidy(statement: str) -> str:
    """SQLAlchemy indents DDL with tabs and leaves trailing spaces on each line.

    Both are harmless, but the rendered files are meant to be read, so they get
    normalised to four spaces with no trailing whitespace.
    """
    lines = [line.replace("\t", "    ").rstrip() for line in statement.strip().splitlines()]
    return "\n".join(lines) + ";"


def analyze_database(session: Session) -> None:
    """Refresh the query planner's table statistics.

    Worth running after a bulk load and not much else. Without stats SQLite
    guesses at join order, and on the 54k row scan it picks the wrong driving
    table for the reference join: the plan falls back to scanning the
    materialised CTE instead of building an automatic covering index over it.
    Measured on that scan, ANALYZE takes the trend query from 0.21s to 0.09s.

    Runs on the session's own connection rather than opening a new one, because
    an in-memory SQLite database is per-connection and a second connection would
    analyze an empty copy.
    """
    session.execute(text("ANALYZE"))
    session.commit()


def render_schema_sql() -> str:
    """DDL for every table, in dependency order.

    Used by scripts/dump_schema_sql.py so sql/schema.sql is always a rendering
    of the models rather than a second copy that can drift.
    """
    dialect = sqlite_dialect.dialect()
    statements = [
        _tidy(str(CreateTable(table).compile(dialect=dialect)))
        for table in Base.metadata.sorted_tables
    ]
    return "\n\n".join(statements) + "\n"


def render_indexes_sql() -> str:
    """DDL for the explicit indexes, grouped by table."""
    dialect = sqlite_dialect.dialect()
    lines: list[str] = []
    for table in Base.metadata.sorted_tables:
        indexes = sorted(table.indexes, key=lambda idx: idx.name or "")
        if not indexes:
            continue
        lines.append(f"-- {table.name}")
        for index in indexes:
            lines.append(_tidy(str(CreateIndex(index).compile(dialect=dialect))))
        lines.append("")
    return "\n".join(lines)
