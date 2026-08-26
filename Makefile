# Thin wrappers over the CLI. Everything here also works as a plain python call,
# the targets just save typing and document the order the stages run in.

PYTHON ?= python
DB ?= sqlite:///data/campaign.db
ENGINE ?= stub
CATALOG ?= data/scan_25nt.csv
TARGET_NAME ?= cdna_target
TARGET_SEQUENCE ?= GCTTCCAGCTTATTGAATTACACGCAGAGGGTAGCGGCTCTGCGCATTCAATTGCTGCGCGCTGAAGCGCGGAAGC
DIMENSIONS ?= truncation paired_bases gc position substitution

.PHONY: help install db import sweep rank strongest trends explain runs demo schema test lint clean

help:
	@echo "make install    install the package and dev extras"
	@echo "make demo       import the scan CSV, rank it, show the trends"
	@echo "make import     load a scan CSV and unpivot its energy columns"
	@echo "make sweep      simulate instead of importing (ENGINE=stub or nupack)"
	@echo "make rank       rank by delta_g_binding, most positive first"
	@echo "make strongest  the mirror ranking, most negative first"
	@echo "make trends     aggregate over DIMENSIONS (see: make trends-list)"
	@echo "make explain    EXPLAIN QUERY PLAN for the ranking query"
	@echo "make runs       list runs, conditions and selections"
	@echo "make schema     regenerate sql/schema.sql and sql/indexes.sql"
	@echo "make test       run the test suite"

install:
	$(PYTHON) -m pip install -e ".[dev]"

db:
	$(PYTHON) -m nupack_campaign.cli --db $(DB) init-db

import: db
	$(PYTHON) -m nupack_campaign.cli --db $(DB) import-results \
		--catalog $(CATALOG) \
		--target-name $(TARGET_NAME) \
		--target-sequence $(TARGET_SEQUENCE)

sweep: db
	$(PYTHON) -m nupack_campaign.cli --db $(DB) load-variants --catalog $(CATALOG)
	$(PYTHON) -m nupack_campaign.cli --db $(DB) sweep \
		--engine $(ENGINE) \
		--target-name $(TARGET_NAME) \
		--target-sequence $(TARGET_SEQUENCE)

rank:
	$(PYTHON) -m nupack_campaign.cli --db $(DB) rank --out data/shortlist.csv

strongest:
	$(PYTHON) -m nupack_campaign.cli --db $(DB) rank --strongest

trends:
	$(PYTHON) -m nupack_campaign.cli --db $(DB) trends --by $(DIMENSIONS) --out data/trends

trends-list:
	$(PYTHON) -m nupack_campaign.cli --db $(DB) trends --list

explain:
	$(PYTHON) -m nupack_campaign.cli --db $(DB) explain

runs:
	$(PYTHON) -m nupack_campaign.cli --db $(DB) runs

demo: clean import rank trends

schema:
	$(PYTHON) scripts/dump_schema_sql.py

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check src tests scripts

clean:
	rm -f data/campaign.db data/campaign.db-wal data/campaign.db-shm data/shortlist.csv
	rm -rf data/trends
