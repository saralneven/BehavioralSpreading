import re
from pathlib import Path
import pandas as pd

# =========================
# USER SETTINGS
# =========================
INPUT_DIR  = Path("../../outputs/simulations/GMM_3D") / "<run_folder_name>"  # update with actual run folder   # folder containing Sim*.csv
OUTPUT_DIR = Path("../../outputs/simulations/GMM_3D/Grid-Analysis")       # where you want merged chunks saved
OUTPUT_FORMAT = "csv.gz"  # "parquet" (recommended) or "csv.gz"

# If True, reads CSVs in chunks to reduce RAM usage (only used for csv.gz output)
CSV_CHUNKSIZE = 250_000


# =========================
# FILENAME PARSER
# =========================
# Example:
# ModelDi_A1p0_N5_NFR1_FRspdgrid_FRv0p400_dec0p0333_20260116-171250_NFR1-2-3-4-5_FRspdgrid_dec0p0333_tmax400
#
# We will extract:
# - model: Di / DiT / DiTDe / DiTDeS etc
# - area:  1p0  (optional, but nice to store)
# - N:     5    (group size in filename, if you encode it there)
# - NFR:   1
#
# IMPORTANT: If your actual saved CSVs have a slightly different prefix/suffix,
# just tweak the regex below.
FNAME_RE = re.compile(
    r"Model(?P<model>[A-Za-z0-9]+)"
    r"_A(?P<area>\d+p\d+)"
    r"_N(?P<N>\d+)"
    r"_NFR(?P<NFR>\d+)"
)

def parse_meta_from_name(fname: str):
    m = FNAME_RE.search(fname)
    if not m:
        return None
    d = m.groupdict()
    d["N"] = int(d["N"])
    d["NFR"] = int(d["NFR"])
    # convert "1p0" -> 1.0
    d["area"] = float(d["area"].replace("p", "."))
    return d


# =========================
# DISCOVER FILES
# =========================
INPUT_DIR = INPUT_DIR.expanduser().resolve()
OUTPUT_DIR = OUTPUT_DIR.expanduser().resolve()
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

csv_files = sorted(INPUT_DIR.rglob("*.csv"))
if not csv_files:
    raise FileNotFoundError(f"No .csv files found under {INPUT_DIR}")

# Group input files by (model, NFR)
groups = {}
skipped = []
for fp in csv_files:
    meta = parse_meta_from_name(fp.stem)
    if meta is None:
        skipped.append(fp.name)
        continue
    key = (meta["model"], meta["NFR"])
    groups.setdefault(key, []).append((fp, meta))

print(f"Found {len(csv_files)} CSVs total.")
print(f"Grouped into {len(groups)} (model, N_FR) chunks.")
if skipped:
    print(f"Skipped {len(skipped)} files (filename didn't match parser). Example: {skipped[0]}")


# =========================
# WRITE MERGED CHUNKS
# =========================
def out_path(model: str, nfr: int):
    sub = OUTPUT_DIR / f"model={model}" / f"NFR={nfr}"
    sub.mkdir(parents=True, exist_ok=True)
    if OUTPUT_FORMAT.lower() == "parquet":
        return sub / "merged.parquet"
    elif OUTPUT_FORMAT.lower() == "csv.gz":
        return sub / "merged.csv.gz"
    else:
        raise ValueError("OUTPUT_FORMAT must be 'parquet' or 'csv.gz'")

for (model, nfr), file_list in groups.items():
    file_list = sorted(file_list, key=lambda x: x[0].name)
    op = out_path(model, nfr)
    print(f"\nMerging: model={model}, NFR={nfr} -> {op}")
    print(f"  Files: {len(file_list)}")

    if OUTPUT_FORMAT.lower() == "parquet":
        # Parquet approach: concatenate in memory per group.
        # If this group is still too big for RAM, tell me and I’ll switch you to
        # a streaming parquet writer (pyarrow dataset) approach.
        dfs = []
        for fp, meta in file_list:
            df = pd.read_csv(fp)
            # add metadata columns from filename (handy for checks)
            df["model"] = meta["model"]
            df["area_label"] = meta["area"]
            df["N_from_fname"] = meta["N"]
            df["NFR_from_fname"] = meta["NFR"]
            df["source_file"] = fp.name
            dfs.append(df)

        merged = pd.concat(dfs, ignore_index=True)
        merged.to_parquet(op, index=False)  # compression is handled by engine defaults; can be tuned if needed
        print(f"  Wrote rows: {len(merged):,}")

    else:
        # csv.gz approach: stream append to avoid holding all data in memory
        # (writes header once, then appends).
        first = True
        total = 0
        for fp, meta in file_list:
            for chunk in pd.read_csv(fp, chunksize=CSV_CHUNKSIZE):
                chunk["model"] = meta["model"]
                chunk["area_label"] = meta["area"]
                chunk["N_from_fname"] = meta["N"]
                chunk["NFR_from_fname"] = meta["NFR"]
                chunk["source_file"] = fp.name

                chunk.to_csv(op, mode="wt" if first else "at",
                             header=first, index=False, compression="gzip")
                first = False
                total += len(chunk)
        print(f"  Wrote rows: {total:,}")

print("\nDone.")
