#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
import socket
import platform
from datetime import datetime, timezone
from influxdb_client_3 import InfluxDBClient3, Point

# Try loading local .env
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

def compute_run_id(json_bytes: bytes) -> str:
    h = hashlib.sha256()
    h.update(json_bytes)
    return h.hexdigest()

def parse_benchmark_name(full_name):
    """
    Splits 'BM_Suite/1024' into ('BM_Suite', '1024').
    Splits 'BM_Simple' into ('BM_Simple', 'N/A').
    """
    if "/" in full_name:
        parts = full_name.split("/", 1)
        return parts[0], parts[1]
    return full_name, "N/A"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True, help="Path to Google Benchmark JSON")
    args = parser.parse_args()

    # 1. Load JSON
    try:
        with open(args.json, "rb") as f:
            js_bytes = f.read()
        data = json.loads(js_bytes.decode("utf-8"))
    except Exception as e:
        print(f"[error] Could not read JSON file: {e}")
        return

    benchmarks = data.get("benchmarks", [])
    if not benchmarks:
        print("[error] JSON file is valid but contains NO benchmarks.")
        return

    context = data.get("context", {})
    run_id = context.get("run_id") or compute_run_id(js_bytes)
    host_name = socket.gethostname()

    # Timestamp logic
    gb_date_str = context.get("date")
    if gb_date_str:
        try:
            # Try parsing standard ISO formats
            timestamp = datetime.strptime(gb_date_str, "%Y-%m-%dT%H:%M:%S%z")
        except:
            # Fallback for simple formats or just use now()
            timestamp = datetime.now(timezone.utc)
    else:
        timestamp = datetime.now(timezone.utc)

    # 2. Connect to InfluxDB
    # Windows Cert Fix
    if platform.system() == "Windows":
        try:
            import certifi
            os.environ["SSL_CERT_FILE"] = certifi.where()
        except ImportError:
            pass

    token = os.environ.get("INFLUXDB_TOKEN")
    host = os.environ.get("INFLUXDB_HOST")
    database = os.environ.get("INFLUXDB_DATABASE")

    if not token or not host or not database:
        print("[error] Missing INFLUXDB environment variables.")
        return

    client_args = {"host": host, "token": token, "database": database}

    points = []

    print(f"--- Parsing {len(benchmarks)} benchmarks ---")

    for b in benchmarks:
        full_name = b.get("name")
        suite, arg = parse_benchmark_name(full_name)

        real_time = b.get("real_time")
        cpu_time = b.get("cpu_time")

        if b.get("error_occurred"):
            continue

        p = (
            Point("ga_benchmark")
            .tag("benchmark", full_name)
            .tag("suite", suite)
            .tag("arg", arg)
            .tag("host", host_name)
            .tag("run_id", run_id)
            .tag("git_branch", context.get("git_branch", "unknown"))
            .time(timestamp)
        )

        # -------------------------------------------------------
        # FIX: Explicit Types to avoid Schema Conflicts
        # -------------------------------------------------------

        # FLOATS
        if real_time is not None: p.field("real_time_ns", float(real_time))
        if cpu_time is not None: p.field("cpu_time_ns", float(cpu_time))

        # INTEGERS (Must match InfluxDB Schema)
        if b.get("iterations") is not None:
            p.field("iterations", int(b.get("iterations")))

        # Memory metrics (Safe to keep as floats usually, but bytes are technically ints)
        # We will use float for averages like allocs_per_iter
        if b.get("allocs_per_iter") is not None:
            p.field("allocs_per_iter", float(b.get("allocs_per_iter")))

        # Bytes are discrete, but keeping as float is safer for averages unless
        # your schema strictly enforces integer for these too.
        # If you get another 400 error for these, change them to int().
        if b.get("max_bytes_used") is not None:
            p.field("max_bytes_used", float(b.get("max_bytes_used")))
        if b.get("net_heap_growth") is not None:
            p.field("net_heap_growth", float(b.get("net_heap_growth")))

        points.append(p)

    # 3. Upload
    if not points:
        print("[upload] No valid points generated.")
        return

    print(f"[upload] Uploading {len(points)} points to InfluxDB...")
    try:
        with InfluxDBClient3(**client_args) as client:
            client.write(record=points)
        print("[upload] Success!")
    except Exception as e:
        print(f"[error] Upload failed: {e}")
        # Optional: Print detailed error if available
        try:
            print(f"Reason: {e.reason}")
            print(f"Body: {e.body}")
        except:
            pass

if __name__ == "__main__":
    main()