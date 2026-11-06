#!/usr/bin/env python3
"""Build a reproducible, entirely synthetic respiratory-virus surveillance DB.

Study draft: all records are simulated. No real people or observations are used.
Run from this folder with: python generate_database.py
The default output is respiratory_surveillance.sqlite beside this script.
"""

from __future__ import annotations

import argparse
import math
import random
import sqlite3
from collections import defaultdict
from datetime import date, timedelta
from pathlib import Path
from typing import Any


DEFAULT_SEED = 20261106
START_DATE = date(2025, 10, 1)
END_DATE = date(2026, 3, 31)
PARTICIPANT_COUNT = 1_500
SPECIMEN_COUNT = 2_400


# Core tables are normalized reference, participant, specimen and assay data.
# A composite FK connects each result to the assay/pathogen pair. Deliberate
# source-feed duplicates live only in raw_lab_feed, where intake rows need a
# surrogate key so questionable source records can be retained for cleaning.
SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE region (
    region_code TEXT PRIMARY KEY,
    region_name TEXT NOT NULL UNIQUE,
    climate_band TEXT NOT NULL CHECK (climate_band IN ('coastal', 'inland', 'upland'))
);

CREATE TABLE site (
    site_id INTEGER PRIMARY KEY,
    region_code TEXT NOT NULL REFERENCES region(region_code),
    site_name TEXT NOT NULL UNIQUE,
    site_type TEXT NOT NULL CHECK (site_type IN ('community clinic', 'outpatient', 'urgent care', 'hospital')),
    catchment_population INTEGER NOT NULL CHECK (catchment_population > 0)
);

CREATE TABLE participant (
    participant_id INTEGER PRIMARY KEY,
    age_band TEXT NOT NULL CHECK (age_band IN ('under 5', '5-17', '18-49', '50-64', '65+')),
    risk_category TEXT NOT NULL CHECK (risk_category IN ('standard', 'higher')),
    enrolment_group TEXT NOT NULL CHECK (enrolment_group IN ('community', 'care setting', 'health workforce'))
);

CREATE TABLE pathogen (
    pathogen_code TEXT PRIMARY KEY,
    pathogen_name TEXT NOT NULL UNIQUE,
    virus_family TEXT NOT NULL,
    severity_order INTEGER NOT NULL CHECK (severity_order BETWEEN 1 AND 4)
);

CREATE TABLE assay (
    assay_code TEXT NOT NULL CHECK (assay_code IN ('PCR', 'ANTIGEN')),
    pathogen_code TEXT NOT NULL REFERENCES pathogen(pathogen_code),
    assay_name TEXT NOT NULL,
    typical_turnaround TEXT NOT NULL CHECK (typical_turnaround IN ('hours', 'minutes')),
    PRIMARY KEY (assay_code, pathogen_code),
    UNIQUE (pathogen_code, assay_name)
);

CREATE TABLE specimen (
    specimen_id INTEGER PRIMARY KEY,
    participant_id INTEGER NOT NULL REFERENCES participant(participant_id),
    site_id INTEGER NOT NULL REFERENCES site(site_id),
    collected_on TEXT NOT NULL,
    symptom_onset TEXT,
    symptomatic INTEGER NOT NULL CHECK (symptomatic IN (0, 1)),
    symptom_severity TEXT NOT NULL CHECK (symptom_severity IN ('none', 'mild', 'moderate', 'severe')),
    temperature_celsius REAL CHECK (temperature_celsius BETWEEN 35.0 AND 41.5),
    CHECK (date(collected_on) IS NOT NULL),
    CHECK (symptom_onset IS NULL OR date(symptom_onset) IS NOT NULL),
    CHECK (symptom_onset IS NULL OR symptom_onset <= collected_on),
    CHECK ((symptomatic = 0 AND symptom_severity = 'none') OR (symptomatic = 1 AND symptom_severity <> 'none'))
);

CREATE TABLE lab_result (
    specimen_id INTEGER NOT NULL REFERENCES specimen(specimen_id),
    assay_code TEXT NOT NULL,
    pathogen_code TEXT NOT NULL,
    replicate_no INTEGER NOT NULL CHECK (replicate_no BETWEEN 1 AND 3),
    result_category TEXT NOT NULL CHECK (result_category IN ('detected', 'not detected', 'inconclusive')),
    ct_value REAL CHECK (ct_value BETWEEN 5.0 AND 45.0),
    viral_load_copies_ml INTEGER CHECK (viral_load_copies_ml >= 0),
    turnaround_hours REAL NOT NULL CHECK (turnaround_hours > 0 AND turnaround_hours <= 72),
    PRIMARY KEY (specimen_id, assay_code, pathogen_code, replicate_no),
    FOREIGN KEY (assay_code, pathogen_code) REFERENCES assay(assay_code, pathogen_code),
    CHECK (result_category = 'detected' OR (ct_value IS NULL AND viral_load_copies_ml IS NULL))
);

CREATE TABLE participant_vaccination (
    participant_id INTEGER NOT NULL REFERENCES participant(participant_id),
    pathogen_code TEXT NOT NULL REFERENCES pathogen(pathogen_code),
    dose_no INTEGER NOT NULL CHECK (dose_no BETWEEN 1 AND 4),
    administered_on TEXT NOT NULL,
    product_category TEXT NOT NULL CHECK (product_category IN ('mRNA', 'protein', 'inactivated', 'seasonal formulation')),
    PRIMARY KEY (participant_id, pathogen_code, dose_no),
    CHECK (date(administered_on) IS NOT NULL)
);

CREATE TABLE weekly_surveillance (
    site_id INTEGER NOT NULL REFERENCES site(site_id),
    pathogen_code TEXT NOT NULL REFERENCES pathogen(pathogen_code),
    week_start TEXT NOT NULL,
    tests_count INTEGER NOT NULL CHECK (tests_count > 0),
    positives_count INTEGER NOT NULL CHECK (positives_count BETWEEN 0 AND tests_count),
    positivity_rate REAL NOT NULL CHECK (positivity_rate BETWEEN 0.0 AND 1.0),
    alert_level TEXT NOT NULL CHECK (alert_level IN ('routine', 'watch', 'elevated')),
    PRIMARY KEY (site_id, pathogen_code, week_start),
    CHECK (date(week_start) IS NOT NULL),
    CHECK (ABS(positivity_rate - ROUND(1.0 * positives_count / tests_count, 5)) <= 0.000001),
    CHECK ((alert_level = 'routine' AND positivity_rate < 0.12)
        OR (alert_level = 'watch' AND positivity_rate >= 0.12 AND positivity_rate < 0.24)
        OR (alert_level = 'elevated' AND positivity_rate >= 0.24))
);

-- Raw intake is deliberately less constrained than curated tables: duplicate
-- external keys and absent receipt timestamps are retained for QC exercises.
CREATE TABLE raw_lab_feed (
    feed_row_id INTEGER PRIMARY KEY,
    source_record_key TEXT NOT NULL,
    specimen_id INTEGER NOT NULL,
    assay_code TEXT NOT NULL,
    pathogen_code TEXT NOT NULL,
    reported_result TEXT NOT NULL,
    source_received_at TEXT,
    ingested_at TEXT NOT NULL
);

CREATE INDEX idx_specimen_date_site ON specimen(collected_on, site_id);
CREATE INDEX idx_result_pathogen_category ON lab_result(pathogen_code, result_category);
CREATE INDEX idx_raw_source_key ON raw_lab_feed(source_record_key);

CREATE VIEW v_outbreak_monitoring AS
SELECT w.week_start, r.region_name, s.site_name, s.site_type,
       p.pathogen_name, w.tests_count, w.positives_count,
       w.positivity_rate, w.alert_level
FROM weekly_surveillance AS w
JOIN site AS s ON s.site_id = w.site_id
JOIN region AS r ON r.region_code = s.region_code
JOIN pathogen AS p ON p.pathogen_code = w.pathogen_code;
"""


REGIONS = [
    ("NORTH", "North District", "upland"),
    ("CENTRAL", "Central District", "inland"),
    ("EAST", "East District", "coastal"),
    ("SOUTH", "South District", "coastal"),
]

SITES = [
    ("North Community Clinic", "NORTH", "community clinic", 28_000),
    ("North General Hospital", "NORTH", "hospital", 82_000),
    ("Central Outpatient Centre", "CENTRAL", "outpatient", 54_000),
    ("Central Urgent Care", "CENTRAL", "urgent care", 37_000),
    ("East Community Clinic", "EAST", "community clinic", 31_000),
    ("East General Hospital", "EAST", "hospital", 91_000),
    ("South Outpatient Centre", "SOUTH", "outpatient", 58_000),
    ("South Urgent Care", "SOUTH", "urgent care", 42_000),
]

PATHOGENS = [
    ("SARS_COV2", "SARS-CoV-2", "Coronaviridae", 3),
    ("FLU_A", "Influenza A", "Orthomyxoviridae", 3),
    ("FLU_B", "Influenza B", "Orthomyxoviridae", 2),
    ("RSV", "Respiratory syncytial virus", "Pneumoviridae", 3),
]

AGE_WEIGHTS = [
    ("under 5", 0.09),
    ("5-17", 0.16),
    ("18-49", 0.39),
    ("50-64", 0.18),
    ("65+", 0.18),
]


def _iso_day(value: date) -> str:
    return value.isoformat()


def _monday(value: date) -> str:
    return _iso_day(value - timedelta(days=value.weekday()))


def _sample_collection_date(rng: random.Random, date_weights: list[float]) -> date:
    offset = rng.choices(range((END_DATE - START_DATE).days + 1), weights=date_weights, k=1)[0]
    return START_DATE + timedelta(days=offset)


def _positivity_probability(pathogen_code: str, collected_on: date, site_type: str) -> float:
    day_index = (collected_on - START_DATE).days
    peak = {"SARS_COV2": 100, "FLU_A": 104, "FLU_B": 132, "RSV": 55}[pathogen_code]
    baseline = {"SARS_COV2": 0.055, "FLU_A": 0.035, "FLU_B": 0.020, "RSV": 0.050}[pathogen_code]
    amplitude = {"SARS_COV2": 0.095, "FLU_A": 0.175, "FLU_B": 0.085, "RSV": 0.145}[pathogen_code]
    seasonal = math.exp(-((day_index - peak) / 53.0) ** 2)
    site_factor = {"community clinic": 0.92, "outpatient": 1.0, "urgent care": 1.12, "hospital": 1.18}[site_type]
    return min(0.48, (baseline + amplitude * seasonal) * site_factor)


def _temperature(rng: random.Random, severity: str) -> float:
    shift = {"mild": 0.15, "moderate": 0.65, "severe": 1.25}[severity]
    return round(min(40.6, max(35.4, rng.gauss(37.0 + shift, 0.55))), 1)


def _draw_result(rng: random.Random, positive_probability: float) -> tuple[str, float | None, int | None]:
    if rng.random() < positive_probability:
        # Log-normal range is deliberately broad. Ct is inversely associated
        # with load, while both remain fictional, course-scale measurements.
        log10_load = min(8.0, max(2.0, rng.gauss(5.0, 1.05)))
        load = int(10**log10_load)
        ct = round(min(42.0, max(12.0, 40.0 - 3.2 * log10_load + rng.gauss(0, 1.2))), 1)
        return "detected", ct, load
    if rng.random() < 0.004:
        return "inconclusive", None, None
    return "not detected", None, None


def _insert_reference_data(connection: sqlite3.Connection, rng: random.Random) -> None:
    connection.executemany("INSERT INTO region VALUES (?, ?, ?)", REGIONS)
    connection.executemany(
        "INSERT INTO site(region_code, site_name, site_type, catchment_population) VALUES (?, ?, ?, ?)",
        [(region, name, kind, pop) for name, region, kind, pop in SITES],
    )
    connection.executemany("INSERT INTO pathogen VALUES (?, ?, ?, ?)", PATHOGENS)
    connection.executemany(
        "INSERT INTO assay VALUES (?, ?, ?, ?)",
        [(assay, code, f"{name} {assay.lower()} assay", "hours" if assay == "PCR" else "minutes")
         for code, name, _family, _severity in PATHOGENS
         for assay in ("PCR", "ANTIGEN")],
    )

    age_values, age_weights = zip(*AGE_WEIGHTS)
    participants = []
    for participant_id in range(1, PARTICIPANT_COUNT + 1):
        age = rng.choices(age_values, weights=age_weights, k=1)[0]
        risk = "higher" if (age in ("under 5", "65+") and rng.random() < 0.60) or rng.random() < 0.15 else "standard"
        enrolment = rng.choices(
            ["community", "care setting", "health workforce"],
            weights=[0.78, 0.09, 0.13], k=1,
        )[0]
        participants.append((participant_id, age, risk, enrolment))
    connection.executemany("INSERT INTO participant VALUES (?, ?, ?, ?)", participants)

    vaccinations = []
    vaccine_start = START_DATE - timedelta(days=250)
    for participant_id, age, _risk, _group in participants:
        for code, _name, _family, _severity in PATHOGENS:
            p_vaccinated = 0.58 if code in ("FLU_A", "FLU_B") else 0.48
            if code == "RSV":
                p_vaccinated = 0.38 if age in ("65+", "50-64") else 0.035
            if rng.random() >= p_vaccinated:
                continue
            doses = 1 if rng.random() < 0.72 else 2
            product = rng.choice(["mRNA", "protein", "inactivated", "seasonal formulation"])
            for dose_no in range(1, doses + 1):
                administered = vaccine_start + timedelta(days=rng.randint(0, (END_DATE - vaccine_start).days))
                vaccinations.append((participant_id, code, dose_no, _iso_day(administered), product))
    connection.executemany("INSERT INTO participant_vaccination VALUES (?, ?, ?, ?, ?)", vaccinations)


def _insert_surveillance(connection: sqlite3.Connection, rng: random.Random) -> None:
    day_count = (END_DATE - START_DATE).days + 1
    date_weights = []
    for offset in range(day_count):
        # Mildly higher collection volume in winter, with weekly noise.
        peak = math.exp(-((offset - 93) / 68.0) ** 2)
        date_weights.append(0.72 + 0.56 * peak + rng.uniform(0.0, 0.18))

    participant_ids = list(range(1, PARTICIPANT_COUNT + 1))
    site_meta = {i + 1: (region, name, kind, pop) for i, (name, region, kind, pop) in enumerate(SITES)}
    participant_meta = {row[0]: row for row in connection.execute(
        "SELECT participant_id, age_band, risk_category, enrolment_group FROM participant"
    )}

    specimens: list[tuple[Any, ...]] = []
    result_rows: list[tuple[Any, ...]] = []
    result_map_for_feed: list[tuple[int, str, str, int, str]] = []
    site_weights = [entry[3] for entry in SITES]

    for specimen_id in range(1, SPECIMEN_COUNT + 1):
        collected = _sample_collection_date(rng, date_weights)
        site_id = rng.choices(range(1, len(SITES) + 1), weights=site_weights, k=1)[0]
        participant_id = rng.choice(participant_ids)
        age, risk, enrolment = participant_meta[participant_id][1:]
        risk_factor = 1.18 if risk == "higher" else 1.0
        age_factor = 1.14 if age in ("under 5", "65+") else 1.0
        symptomatic = 1 if rng.random() < 0.83 else 0
        severity = "none"
        onset: str | None = None
        temp: float | None = None
        if symptomatic:
            severity = rng.choices(["mild", "moderate", "severe"], weights=[0.67, 0.27, 0.06], k=1)[0]
            # About 15% of onset dates are intentionally absent in this curated
            # table; among the remainder onset is 0-7 days before collection.
            if rng.random() >= 0.15:
                onset = _iso_day(collected - timedelta(days=rng.randint(0, 7)))
            # Temperature is not always recorded at intake (intentional NULLs).
            if rng.random() >= 0.06:
                temp = _temperature(rng, severity)
        else:
            # Asymptomatic records can still lack an onset date; fever measure
            # is not applicable and is represented as NULL rather than zero.
            if rng.random() >= 0.92:
                temp = round(rng.uniform(36.2, 37.3), 1)

        specimens.append((
            specimen_id, participant_id, site_id, _iso_day(collected), onset,
            symptomatic, severity, temp,
        ))

        site_type = site_meta[site_id][2]
        for pathogen_code, _name, _family, _severity_order in PATHOGENS:
            probability = _positivity_probability(pathogen_code, collected, site_type) * risk_factor * age_factor
            pcr_result, ct, load = _draw_result(rng, min(probability, 0.62))
            # A small number of positives have intentionally unrecorded load or
            # Ct values. For negatives and inconclusives, both are not applicable.
            if pcr_result == "detected" and rng.random() < 0.055:
                load = None
            if pcr_result == "detected" and rng.random() < 0.035:
                ct = None
            pcr_hours = round(min(36.0, max(2.0, rng.lognormvariate(math.log(10.0), 0.45))), 2)
            result_rows.append((specimen_id, "PCR", pathogen_code, 1, pcr_result, ct, load, pcr_hours))
            result_map_for_feed.append((specimen_id, "PCR", pathogen_code, 1, pcr_result))

            # Repeats and antigen confirmation make laboratory activity uneven,
            # as it would be in a surveillance workflow.
            if rng.random() < 0.035:
                repeat_ct = ct if pcr_result == "detected" else None
                repeat_load = load if pcr_result == "detected" else None
                result_rows.append((specimen_id, "PCR", pathogen_code, 2, pcr_result, repeat_ct, repeat_load,
                                    round(rng.uniform(3.0, 24.0), 2)))
                result_map_for_feed.append((specimen_id, "PCR", pathogen_code, 2, pcr_result))
            if pcr_result == "detected" and rng.random() < 0.22:
                antigen_result = "detected" if rng.random() < 0.86 else "not detected"
                antigen_ct = None
                antigen_load = None
                antigen_hours = round(rng.uniform(0.2, 2.5), 2)
                result_rows.append((specimen_id, "ANTIGEN", pathogen_code, 1, antigen_result,
                                    antigen_ct, antigen_load, antigen_hours))
                result_map_for_feed.append((specimen_id, "ANTIGEN", pathogen_code, 1, antigen_result))

    connection.executemany("INSERT INTO specimen VALUES (?, ?, ?, ?, ?, ?, ?, ?)", specimens)
    connection.executemany("INSERT INTO lab_result VALUES (?, ?, ?, ?, ?, ?, ?, ?)", result_rows)

    # This intentionally duplicated source feed demonstrates a realistic
    # deduplication task without weakening keys on curated analytical tables.
    raw_rows = []
    for feed_id, (specimen_id, assay_code, pathogen_code, replicate_no, result) in enumerate(result_map_for_feed, start=1):
        source_key = f"EXT-{specimen_id:05d}-{pathogen_code}-{assay_code}-{replicate_no}"
        received_at = None if rng.random() < 0.035 else f"{specimens[specimen_id - 1][3]} {rng.randint(6, 22):02d}:{rng.randint(0, 59):02d}:00"
        row = (feed_id, source_key, specimen_id, assay_code, pathogen_code, result, received_at,
               f"{specimens[specimen_id - 1][3]} 23:30:00")
        raw_rows.append(row)
        # Insert an exact repeated source event for a seeded 2.5% of feed rows.
        if rng.random() < 0.025:
            feed_id += 1
            raw_rows.append((feed_id, *row[1:]))
    # Re-number the staging surrogate key after duplicates so it stays unique.
    raw_rows = [(i + 1, *row[1:]) for i, row in enumerate(raw_rows)]
    connection.executemany("INSERT INTO raw_lab_feed VALUES (?, ?, ?, ?, ?, ?, ?, ?)", raw_rows)

    # A materialized weekly roll-up supports fast alert review. It is derived
    # exclusively from lab_result, so the validation file can reconcile it.
    aggregates: dict[tuple[int, str, str], list[int]] = defaultdict(lambda: [0, 0])
    for site_id, collected_on, pathogen_code, result in connection.execute(
        """SELECT s.site_id, s.collected_on, lr.pathogen_code, lr.result_category
           FROM specimen AS s JOIN lab_result AS lr ON lr.specimen_id = s.specimen_id"""
    ):
        key = (site_id, pathogen_code, _monday(date.fromisoformat(collected_on)))
        aggregates[key][0] += 1
        if result == "detected":
            aggregates[key][1] += 1

    summary_rows = []
    for (site_id, pathogen_code, week_start), (tests, positives) in sorted(aggregates.items()):
        rate = positives / tests if tests else 0.0
        stored_rate = round(rate, 5)
        alert = "elevated" if stored_rate >= 0.24 else "watch" if stored_rate >= 0.12 else "routine"
        summary_rows.append((site_id, pathogen_code, week_start, tests, positives, stored_rate, alert))
    connection.executemany("INSERT INTO weekly_surveillance VALUES (?, ?, ?, ?, ?, ?, ?)", summary_rows)


def build_database(output_path: str | Path, seed: int = DEFAULT_SEED) -> dict[str, Any]:
    """Create the database at output_path and return counts for the README/report."""
    output = Path(output_path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()
    rng = random.Random(seed)
    connection = sqlite3.connect(output)
    try:
        connection.execute("PRAGMA foreign_keys = ON")
        connection.executescript(SCHEMA_SQL)
        with connection:
            _insert_reference_data(connection, rng)
            _insert_surveillance(connection, rng)
        connection.commit()

        fk_issues = list(connection.execute("PRAGMA foreign_key_check"))
        if fk_issues:
            raise RuntimeError(f"Foreign-key validation failed: {fk_issues[:5]}")
        quick_check = connection.execute("PRAGMA quick_check").fetchone()[0]
        if quick_check != "ok":
            raise RuntimeError(f"SQLite quick_check failed: {quick_check}")

        table_names = [row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
        )]
        counts = {name: connection.execute(f'SELECT COUNT(*) FROM "{name}"').fetchone()[0]
                  for name in table_names}
        if counts.get("specimen", 0) < 1_000:
            raise RuntimeError("The specimen table must contain at least 1,000 rows.")

        quality = {
            "specimens_missing_symptom_onset": connection.execute(
                "SELECT COUNT(*) FROM specimen WHERE symptom_onset IS NULL"
            ).fetchone()[0],
            "symptomatic_missing_symptom_onset": connection.execute(
                "SELECT COUNT(*) FROM specimen WHERE symptomatic=1 AND symptom_onset IS NULL"
            ).fetchone()[0],
            "asymptomatic_onset_not_applicable": connection.execute(
                "SELECT COUNT(*) FROM specimen WHERE symptomatic=0 AND symptom_onset IS NULL"
            ).fetchone()[0],
            "specimens_missing_temperature": connection.execute(
                "SELECT COUNT(*) FROM specimen WHERE temperature_celsius IS NULL"
            ).fetchone()[0],
            "symptomatic_missing_temperature": connection.execute(
                "SELECT COUNT(*) FROM specimen WHERE symptomatic=1 AND temperature_celsius IS NULL"
            ).fetchone()[0],
            "asymptomatic_temperature_not_applicable": connection.execute(
                "SELECT COUNT(*) FROM specimen WHERE symptomatic=0 AND temperature_celsius IS NULL"
            ).fetchone()[0],
            "positive_pcr_missing_ct": connection.execute(
                "SELECT COUNT(*) FROM lab_result WHERE assay_code='PCR' AND result_category='detected' AND ct_value IS NULL"
            ).fetchone()[0],
            "positive_pcr_missing_load": connection.execute(
                "SELECT COUNT(*) FROM lab_result WHERE assay_code='PCR' AND result_category='detected' AND viral_load_copies_ml IS NULL"
            ).fetchone()[0],
            "raw_feed_duplicate_key_groups": connection.execute(
                "SELECT COUNT(*) FROM (SELECT source_record_key FROM raw_lab_feed GROUP BY source_record_key HAVING COUNT(*) > 1)"
            ).fetchone()[0],
            "raw_feed_excess_duplicate_rows": connection.execute(
                "SELECT COALESCE(SUM(n-1),0) FROM (SELECT COUNT(*) AS n FROM raw_lab_feed GROUP BY source_record_key HAVING n > 1)"
            ).fetchone()[0],
            "raw_feed_missing_received_at": connection.execute(
                "SELECT COUNT(*) FROM raw_lab_feed WHERE source_received_at IS NULL"
            ).fetchone()[0],
        }
        return {"seed": seed, "path": str(output), "counts": counts, "quality": quality,
                "quick_check": quick_check, "foreign_key_issues": len(fk_issues)}
    finally:
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("respiratory_surveillance.sqlite"),
                        help="destination SQLite file (replaced if it already exists)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help=f"random seed (default: {DEFAULT_SEED})")
    args = parser.parse_args()
    summary = build_database(args.output, args.seed)
    print(f"Created synthetic SQLite database: {summary['path']}")
    print(f"Seed: {summary['seed']}; SQLite quick_check: {summary['quick_check']}; FK issues: {summary['foreign_key_issues']}")
    for table, count in summary["counts"].items():
        print(f"  {table:24s} {count:>6,}")
    print("Intentional quality issues:")
    for issue, count in summary["quality"].items():
        print(f"  {issue:38s} {count:>6,}")


if __name__ == "__main__":
    main()
