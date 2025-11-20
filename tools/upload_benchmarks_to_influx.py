#!/usr/bin/env python3
"""
Upload Google Benchmark JSON results into InfluxDB 3.x using influxdb3-python.

Requires:
    pip install influxdb3-python

Environment variables:
    INFLUXDB_HOST      (e.g. "https://eu-central-1-1.aws.cloud2.influxdata.com")
    INFLUXDB_TOKEN
    INFLUXDB_ORG
    INFLUXDB_DATABASE   (bucket name on cloud, or DB name for local v3)

Usage:
    upload_benchmarks_to_influx.py --json path/to/file.json
"""

import argparse
import hashlib
import json
import os
import socket
from datetime import datetime, timezone

from influxdb_client_3 import InfluxDBClient3, Point


UPLOADED_RUNS_FILE = os.path.expanduser("~/.ga_smith_uploaded_runs")


def load_uploaded_run_ids():
    try:
        with open(UPLOADED_RUNS_FILE, "r", encoding="utf-8") as f:
            return set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        return set()


def save_uploaded_run_ids(run_ids):
    os.makedirs(os.path.dirname(UPLOADED_RUNS_FILE), exist_ok=True)
    with open(UPLOADED_RUNS_FILE, "w", encoding="utf-8") as f:
        for rid in sorted(run_ids):
            f.write(rid + "\n")


def compute_run_id(json_bytes: bytes) -> str:
    h = hashlib.sha256()
    h.update(json_bytes)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True, help="Path to Google Benchmark JSON")
    args = parser.parse_args()

    # -----------------------
    # Load JSON
    # -----------------------
    with open(args.json, "rb") as f:
        js_bytes = f.read()

    run_id = compute_run_id(js_bytes)
    uploaded = load_uploaded_run_ids()

    if run_id in uploaded:
        print(f"[upload] Already uploaded run_id={run_id[:8]}, skipping.")
        return

    data = json.loads(js_bytes.decode("utf-8"))
    context = data.get("context", {})
    benchmarks = data.get("benchmarks", [])

    if not benchmarks:
        print("[upload] No benchmark entries found in JSON.")
        return

    # -----------------------
    # Benchmark metadata
    # -----------------------
    host = socket.gethostname()
    cpu_mhz = context.get("mhz_per_cpu")
    num_cpus = context.get("num_cpus")
    run_id_ctx = context.get("run_id") or run_id

    # Custom context injected by C++
    build_type = context.get("build_type", "unknown")
    compiler = context.get("compiler", "unknown")
    ga_signature = context.get("ga_signature", "unknown")
    git_sha = context.get("git_sha", "unknown")
    git_branch = context.get("git_branch", "unknown")

    # Google Benchmark timestamp
    gb_date_str = context.get("date")
    if gb_date_str:
        try:
            timestamp = datetime.strptime(gb_date_str, "%Y/%m/%d-%H:%M:%S")
            timestamp = timestamp.replace(tzinfo=timezone.utc)
        except Exception:
            timestamp = datetime.now(timezone.utc)
    else:
        timestamp = datetime.now(timezone.utc)

    # -----------------------
    # Connect to InfluxDB 3.x
    # -----------------------
    client = InfluxDBClient3(
        host=os.environ.get("INFLUXDB_HOST"),
        token=os.environ.get("INFLUXDB_TOKEN"),
        org=os.environ.get("INFLUXDB_ORG"),
        database=os.environ.get("INFLUXDB_DATABASE"),
    )

    points = []

    # -----------------------
    # Build points
    # -----------------------
    for b in benchmarks:
        name = b.get("name")
        real_time = b.get("real_time")
        cpu_time = b.get("cpu_time")
        iterations = b.get("iterations")

        # Memory fields from MemoryManager
        allocs_per_iter = b.get("allocs_per_iter")
        max_bytes_used = b.get("max_bytes_used")
        total_allocated_bytes = b.get("total_allocated_bytes")
        net_heap_growth = b.get("net_heap_growth")

        # Simple machine normalization
        cycles_estimate = None
        normalized_time_1ghz_ns = None

        if real_time is not None and cpu_mhz is not None:
            cycles_estimate = float(real_time) * float(cpu_mhz) * 1e-3
            normalized_time_1ghz_ns = cycles_estimate

        p = (
            Point("ga_benchmark")
            .tag("benchmark", name)
            .tag("host", host)
            .tag("run_id", run_id_ctx)
            .tag("build_type", build_type)
            .tag("compiler", compiler)
            .tag("ga_signature", ga_signature)
            .tag("git_sha", git_sha)
            .tag("git_branch", git_branch)
            .time(timestamp)
        )

        if cpu_mhz is not None:
            p.tag("cpu_mhz", str(cpu_mhz))
        if num_cpus is not None:
            p.tag("num_cpus", str(num_cpus))

        # Standard timing fields
        if real_time is not None:
            p.field("real_time_ns", float(real_time))
        if cpu_time is not None:
            p.field("cpu_time_ns", float(cpu_time))
        if iterations is not None:
            p.field("iterations", int(iterations))

        # Memory fields
        if allocs_per_iter is not None:
            p.field("allocs_per_iter", float(allocs_per_iter))
        if max_bytes_used is not None:
            p.field("max_bytes_used", float(max_bytes_used))
        if total_allocated_bytes is not None:
            p.field("total_allocated_bytes", float(total_allocated_bytes))
        if net_heap_growth is not None:
            p.field("net_heap_growth", float(net_heap_growth))

        # Normalized metrics
        if cycles_estimate is not None:
            p.field("cycles_estimate", float(cycles_estimate))
        if normalized_time_1ghz_ns is not None:
            p.field("normalized_time_1ghz_ns", float(normalized_time_1ghz_ns))

        points.append(p)

    # -----------------------
    # Write to InfluxDB 3.x
    # -----------------------
    client.write(points=points)
    client.close()

    # Mark as uploaded
    uploaded.add(run_id)
    save_uploaded_run_ids(uploaded)

    print(f"[upload] Uploaded {len(points)} points (run_id={run_id[:8]}...).")


if __name__ == "__main__":
    main()
