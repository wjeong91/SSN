#!/usr/bin/env python3
"""Download form II Rubisco structures into rubisco_form2/structures/.

- PDB: IDs in ids/pdb.txt plus every PDB entry whose sequence is >= 45% identical to
  9I27 chain A (RCSB sequence search), as gzipped mmCIF
- AlphaFold DB: UniProt accessions in ids/afdb.txt
- ESM Atlas: MGYP IDs in ids/esm.txt

Standard library only, so it runs on a bare GitHub Actions runner.
"""
import gzip
import json
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
IDS = HERE / "ids"
OUT = HERE / "structures"
IDENTITY = 0.45
UA = {"User-Agent": "rubisco-form2-fetch/1.0"}


def get(url, data=None, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={**UA, **({"Content-Type": "application/json"} if data else {})})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001 - retry any network error
            if i == tries - 1:
                raise
            print(f"  retry {url}: {e}")
            time.sleep(2 ** (i + 1))


def read_ids(name):
    return [x.strip() for x in (IDS / name).read_text().split() if x.strip()]


def rcsb_form2_entries():
    seq = "".join(l.strip() for l in (IDS / "9I27_A.fasta").read_text().splitlines() if not l.startswith(">"))
    query = {
        "query": {"type": "terminal", "service": "sequence",
                  "parameters": {"evalue_cutoff": 1e-10, "identity_cutoff": IDENTITY,
                                 "sequence_type": "protein", "value": seq}},
        "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True, "results_verbosity": "verbose"},
    }
    raw = get("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(query).encode())
    (OUT / "rcsb_sequence_search.json").write_bytes(raw)
    hits = json.loads(raw).get("result_set", [])
    return sorted({h["identifier"].split("_")[0].upper() for h in hits})


def save_gz(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        fh.write(data)


def main():
    OUT.mkdir(exist_ok=True)
    manifest = []

    pdb = set(read_ids("pdb.txt"))
    try:
        found = rcsb_form2_entries()
        print(f"RCSB sequence search (>= {IDENTITY:.0%} identity to 9I27_A): {len(found)} entries")
        pdb |= set(found)
    except Exception as e:  # noqa: BLE001
        print(f"RCSB search failed, using listed IDs only: {e}")
    for pid in sorted(pdb):
        path = OUT / "pdb" / f"{pid}.cif.gz"
        try:
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(get(f"https://files.rcsb.org/download/{pid}.cif.gz"))
            manifest.append(("pdb", pid, str(path.relative_to(HERE)), "ok"))
        except Exception as e:  # noqa: BLE001
            manifest.append(("pdb", pid, "", f"failed: {e}"))

    for acc in read_ids("afdb.txt"):
        path = OUT / "afdb" / f"{acc}.cif.gz"
        try:
            if not path.exists():
                entry = json.loads(get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"))[0]
                save_gz(path, get(entry["cifUrl"]))
            manifest.append(("afdb", acc, str(path.relative_to(HERE)), "ok"))
        except Exception as e:  # noqa: BLE001
            manifest.append(("afdb", acc, "", f"failed: {e}"))

    for mid in read_ids("esm.txt"):
        path = OUT / "esm" / f"{mid}.pdb.gz"
        try:
            if not path.exists():
                save_gz(path, get(f"https://api.esmatlas.com/fetchPredictedStructure/{mid}.pdb"))
            manifest.append(("esm", mid, str(path.relative_to(HERE)), "ok"))
        except Exception as e:  # noqa: BLE001
            manifest.append(("esm", mid, "", f"failed: {e}"))

    with open(OUT / "manifest.tsv", "w") as fh:
        fh.write("source\tid\tfile\tstatus\n")
        for row in manifest:
            fh.write("\t".join(row) + "\n")
    ok = sum(r[3] == "ok" for r in manifest)
    print(f"{ok}/{len(manifest)} structures downloaded")
    for r in manifest:
        if r[3] != "ok":
            print("  ", *r)


if __name__ == "__main__":
    main()
