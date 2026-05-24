import argparse
import csv
import gzip
import os
import json
import logging
import re
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, Iterator, Optional, Set, Tuple
import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger("mitomap_cmp")

# Clasificación clínica MitoMap
PATHOGENIC = {"cfrm-p", "cfrm-lp"}
# De menor a mayor prioridad clínica
CLINICAL_PRIORITY = [
    "unclassified",
    "unclear",
    "reported-b",
    "reported-lb",
    "reported",
    "reported-vus",
    "cfrm-vus",
    "conflicting",
    "reported-lp",
    "cfrm-lp",
    "reported-p",
    "cfrm-p",
]

def norm_diseasestatus(x: str) -> Set[str]:
    """
    Normaliza DiseaseStatus a un conjunto de etiquetas canónicas.
    Maneja múltiples valores separados por '|', ',', ';'.
    Ejemplos:
        'Cfrm-[P]'                      -> {'cfrm-p'}
        'Conflicting-reports'           -> {'conflicting'}
        'Reported-/-Unclear'            -> {'reported', 'unclear'}
        'Reported-possibly-synergistic' -> {'reported'}
        '.'                             -> set()
    """
    if not x or x in (".", "None", "nan"):
        return set()
    parts = re.split(r"[|,;]", x)
    out: Set[str] = set()

    for p in parts:
        p = p.strip()
        if not p:
            continue
        # Separar subvalores unidos por ' / ' o '-/-'
        subparts = re.split(r"\s*/\s*|-/-", p)
        for sp in subparts:
            sp = sp.strip().lower()
            # Eliminar corchetes para normalizar primero
            sp_clean = re.sub(r"[\[\]]", "", sp)
            if re.search(r"conflicting", sp_clean):
                out.add("conflicting")
            elif re.search(r"cfrm.*\bp\b", sp_clean) and not re.search(r"lp", sp_clean):
                out.add("cfrm-p")
            elif re.search(r"cfrm.*lp", sp_clean):
                out.add("cfrm-lp")
            elif re.search(r"cfrm.*vus", sp_clean):
                out.add("cfrm-vus")
            elif re.search(r"reported.*\blp\b", sp_clean):
                out.add("reported-lp")
            elif re.search(r"reported.*\bp\b", sp_clean) and not re.search(r"lp|population|possibly|protective", sp_clean):
                out.add("reported-p")
            elif re.search(r"reported.*\bvus\b|\bvus\*\b", sp_clean):
                out.add("reported-vus")
            elif re.search(r"reported.*\blb\b", sp_clean):
                out.add("reported-lb")
            elif re.search(r"reported.*\bb\b", sp_clean):
                out.add("reported-b")
            elif re.search(r"reported", sp_clean):
                out.add("reported")
            elif re.search(r"unclear", sp_clean):
                out.add("unclear")
            else:
                out.add(sp_clean)
    return out

def classify_set(s: Set[str]) -> str:
    """Devuelve la clasificación de mayor prioridad del conjunto."""
    if not s:
        return "unclassified"
    best = "unclassified"
    for v in s:
        if v in CLINICAL_PRIORITY and CLINICAL_PRIORITY.index(v) > CLINICAL_PRIORITY.index(best):
            best = v
    return best

def is_pathogenic(s: Set[str]) -> bool:
    return bool(s & PATHOGENIC)

def is_only_conflicting(s: Set[str]) -> bool:
    return s == {"conflicting"}

def is_classified(x: str) -> bool:
    return bool(x) and x not in (".", "None", "nan", "")

# Extracción de campos VCF MitoMap
def _get_scalar(rec, key: str) -> str:
    v = rec.info.get(key)
    if v is None:
        return "."
    if isinstance(v, (tuple, list)):
        return "|".join(str(i) for i in v if i is not None)
    return str(v)

def get_disease(rec) -> str:
    return _get_scalar(rec, "Disease")

def get_diseasestatus(rec) -> str:
    return _get_scalar(rec, "DiseaseStatus")

def get_aachange(rec) -> str:
    return _get_scalar(rec, "aachange")

def get_heteroplasmy(rec) -> str:
    return _get_scalar(rec, "heteroplasmy")

def get_homoplasmy(rec) -> str:
    return _get_scalar(rec, "homoplasmy")

def get_hgfl(rec) -> str:
    return _get_scalar(rec, "HGFL")

def get_pubmed(rec) -> str:
    return _get_scalar(rec, "PubmedIDs")

def _fmt_gt(gt) -> str:
    if gt is None:
        return "."
    alleles = []
    for a in gt:
        alleles.append("." if a is None else str(a))
    return "/".join(alleles)

def vcf_to_tsv(vcf_file: str, out_file: str, sort_tmp_dir: str = None):
    log.info(f"VCF → TSV: {vcf_file}")
    vcf = pysam.VariantFile(vcf_file)
    tmp_dir = sort_tmp_dir or tempfile.gettempdir()
    tmp_body = tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=tmp_dir, delete=False
    )
    tmp_body_path = tmp_body.name
    header_written = False

    try:
        writer = None
        for rec in vcf.fetch():
            row = {
                "KEY":          f"{rec.chrom}:{rec.pos}:{rec.ref}:{','.join(rec.alts or ['.'])}",
                "CHROM":        rec.chrom,
                "POS":          rec.pos,
                "REF":          rec.ref,
                "ALT":          ",".join(rec.alts or ["."]),
                "DISEASE":      get_disease(rec),
                "DISEASESTATUS": get_diseasestatus(rec),
                "AACHANGE":     get_aachange(rec),
                "HETEROPLASMY": get_heteroplasmy(rec),
                "HOMOPLASMY":   get_homoplasmy(rec),
                "HGFL":         get_hgfl(rec),
                "PUBMED":       get_pubmed(rec),
            }

            for s in rec.samples:
                row[f"{s}_GT"] = _fmt_gt(rec.samples[s].get("GT", None))
            if writer is None:
                writer = csv.DictWriter(
                    tmp_body, fieldnames=list(row.keys()),
                    delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                header_written = True
            writer.writerow(row)
        tmp_body.close()
        if not header_written:
            log.warning(f"VCF vacío: {vcf_file}")
            with gzip.open(out_file, "wt") as f:
                f.write("")
            return

        tmp_sorted_path = tmp_body_path + ".sorted"
        sort_env = os.environ.copy()
        if sort_tmp_dir:
            sort_env["TMPDIR"] = sort_tmp_dir

        sort_cmd = (
            f"(head -n1 {tmp_body_path} && tail -n+2 {tmp_body_path} "
            f"| sort --stable -t$'\\t' -k1,1 "
            f"{'-T ' + sort_tmp_dir if sort_tmp_dir else ''}) "
            f"> {tmp_sorted_path}"
        )
        log.info("Ordenando con sort del sistema...")
        ret = subprocess.run(sort_cmd, shell=True, env=sort_env)
        if ret.returncode != 0:
            raise RuntimeError(f"sort falló con código {ret.returncode}")

        log.info(f"Comprimiendo a {out_file}...")
        with open(tmp_sorted_path, "rb") as src, gzip.open(out_file, "wb") as dst:
            for chunk in iter(lambda: src.read(1 << 20), b""):
                dst.write(chunk)

    finally:
        for p in (tmp_body_path, tmp_body_path + ".sorted"):
            try:
                os.unlink(p)
            except FileNotFoundError:
                pass

def iter_sorted_tsv(file: str) -> Iterator[Tuple[str, dict]]:
    opener = gzip.open if file.endswith(".gz") else open
    with opener(file, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row["KEY"], row

def _sample_cols(row: dict) -> Dict[str, str]:
    return {k[:-3]: v for k, v in row.items() if k.endswith("_GT")}

# Procesado de variante común
def _process_common(
    key: str,
    o: dict,
    n: dict,
    transitions: Dict[Tuple[str, str], int],
    changes_by_transition: Dict[Tuple[str, str], list],
    changed_rows: list,
    counters: dict,
):
    o_raw = o.get("DISEASESTATUS", ".")
    n_raw = n.get("DISEASESTATUS", ".")
    o_set = norm_diseasestatus(o_raw)
    n_set = norm_diseasestatus(n_raw)
    o_label = classify_set(o_set)
    n_label = classify_set(n_set)
    transitions[(o_label, n_label)] += 1

    if o_set != n_set:
        if is_pathogenic(n_set) and not is_pathogenic(o_set) and not is_only_conflicting(n_set):
            sample_data = _sample_cols(n)
            changed_rows.append((key, o_label, n_label,
                                  n.get("DISEASE", "."),
                                  n.get("AACHANGE", "."),
                                  n.get("HETEROPLASMY", "."),
                                  n.get("HOMOPLASMY", "."),
                                  n.get("PUBMED", "."),
                                  sample_data))
            changes_by_transition[(o_label, n_label)].append(key)

        if is_pathogenic(n_set) and not is_pathogenic(o_set) and not is_only_conflicting(n_set):
            counters["gain_path"] += 1
            if not is_classified(o_raw):
                counters["newly_annotated"] += 1
            else:
                counters["reclassified"] += 1

        if is_pathogenic(o_set) and not is_pathogenic(n_set):
            counters["loss_path"] += 1


def compare(old_file: str, new_file: str, outdir: str, sid: str):
    os.makedirs(outdir, exist_ok=True)

    log.info("Comparando con merge externo (streaming)...")

    old_iter = iter_sorted_tsv(old_file)
    new_iter = iter_sorted_tsv(new_file)

    old_current: Optional[Tuple[str, dict]] = next(old_iter, None)
    new_current: Optional[Tuple[str, dict]] = next(new_iter, None)

    n_old = n_new = 0
    only_old = only_new = 0
    common = 0
    old_classified = new_classified = 0

    transitions: Dict[Tuple[str, str], int] = defaultdict(int)
    changes_by_transition: Dict[Tuple[str, str], list] = defaultdict(list)
    changed_rows = []
    counters = dict(gain_path=0, loss_path=0, newly_annotated=0, reclassified=0)
    sample_col_names: list = []

    while old_current is not None or new_current is not None:
        ok = old_current[0] if old_current else None
        nk = new_current[0] if new_current else None

        if ok == nk:
            ov = old_current[1]
            nv = new_current[1]
            n_old += 1; n_new += 1; common += 1

            if is_classified(ov.get("DISEASESTATUS", ".")):
                old_classified += 1
            if is_classified(nv.get("DISEASESTATUS", ".")):
                new_classified += 1

            _process_common(ok, ov, nv, transitions, changes_by_transition,
                            changed_rows, counters)

            if not sample_col_names and changed_rows:
                sample_col_names = list(changed_rows[-1][-1].keys())

            old_current = next(old_iter, None)
            new_current = next(new_iter, None)

        elif nk is None or (ok is not None and ok < nk):
            ov = old_current[1]
            n_old += 1; only_old += 1
            if is_classified(ov.get("DISEASESTATUS", ".")):
                old_classified += 1
            old_current = next(old_iter, None)

        else:
            nv = new_current[1]
            n_new += 1; only_new += 1
            if is_classified(nv.get("DISEASESTATUS", ".")):
                new_classified += 1
            new_current = next(new_iter, None)

    log_file     = os.path.join(outdir, f"{sid}.log")
    json_file    = os.path.join(outdir, f"{sid}.json")
    changed_file = os.path.join(outdir, f"{sid}.changed.tsv")

    with open(changed_file, "w") as f:
        header_cols = ["VARIANT", "OLD", "NEW", "DISEASE",
                       "AACHANGE", "HETEROPLASMY", "HOMOPLASMY", "PUBMED"]
        header_cols.extend(sample_col_names)
        f.write("\t".join(header_cols) + "\n")
        for entry in changed_rows:
            key, o_lbl, n_lbl, disease, aachange, het, hom, pubmed, sample_data = entry
            row_cols = [key, o_lbl, n_lbl, disease, aachange, het, hom, pubmed]
            for col in sample_col_names:
                row_cols.append(str(sample_data.get(col, ".")))
            f.write("\t".join(row_cols) + "\n")

    with open(log_file, "w") as f:
        f.write("RESUMEN DE COMPARACION DE ANOTACIONES MITOMAP\n\n")
        f.write(f"Variantes en OLD: {n_old}\n")
        f.write(f"Variantes en NEW: {n_new}\n\n")
        f.write(f"Variantes en ambas: {common}\n")
        f.write(f"Variantes nuevas (solo NEW): {only_new}\n")
        f.write(f"Variantes perdidas (solo OLD): {only_old}\n\n")

        f.write("Cobertura DiseaseStatus:\n")
        if n_old:
            f.write(f"  OLD clasificados: {old_classified}/{n_old} "
                    f"({old_classified/n_old*100:.2f}%)\n")
        if n_new:
            f.write(f"  NEW clasificados: {new_classified}/{n_new} "
                    f"({new_classified/n_new*100:.2f}%)\n\n")

        f.write(f"Ganancia de patogenicidad (-> Cfrm-P / Cfrm-LP): {counters['gain_path']}\n")
        f.write(f"  Sin clasificacion previa: {counters['newly_annotated']}\n")
        f.write(f"  Reclasificadas: {counters['reclassified']}\n")
        f.write(f"Perdida de patogenicidad: {counters['loss_path']}\n")

        f.write("\n Cambios detallados \n")
        real_changes = {k: v for k, v in changes_by_transition.items() if k[0] != k[1]}
        if not real_changes:
            f.write("  Sin cambios de clasificacion.\n")
        else:
            for (a, b), variants in sorted(real_changes.items(),
                                           key=lambda x: -len(x[1])):
                f.write(f"\n  {a} -> {b}  ({len(variants)} variantes)\n")
                for vk in variants:
                    f.write(f"    {vk}\n")

        cats_present = sorted(
            {a for (a, b) in transitions} | {b for (a, b) in transitions},
            key=lambda x: CLINICAL_PRIORITY.index(x) if x in CLINICAL_PRIORITY else -1
        )
        if cats_present:
            col_w = max(len(c) for c in cats_present) + 2
            row_label_w = col_w
            f.write("\n Matriz de cambios (Filas=OLD, Columnas=NEW)\n\n")
            header_str = " " * (row_label_w + 2)
            for c in cats_present:
                header_str += c.rjust(col_w)
            f.write(f"  {header_str}\n")
            f.write(f"  {' ' * (row_label_w + 2)}{'-' * (col_w * len(cats_present))}\n")
            for a in cats_present:
                row_str = a.ljust(row_label_w) + " |"
                for b in cats_present:
                    row_str += str(transitions.get((a, b), 0)).rjust(col_w)
                f.write(f"  {row_str}\n")

    with open(json_file, "w") as f:
        json.dump({
            "total_old": n_old,
            "total_new": n_new,
            "only_old":  only_old,
            "only_new":  only_new,
            "common":    common,
            "diseasestatus_coverage": {
                "old_classified": old_classified,
                "old_total":      n_old,
                "new_classified": new_classified,
                "new_total":      n_new,
            },
            "gain_pathogenic":            counters["gain_path"],
            "newly_annotated_pathogenic": counters["newly_annotated"],
            "reclassified_pathogenic":    counters["reclassified"],
            "loss_pathogenic":            counters["loss_path"],
            "transitions": {
                f"{a}->{b}": {"count": c, "variants": changes_by_transition.get((a, b), [])}
                for (a, b), c in transitions.items()
            },
            "changed_variants": [
                {
                    "variant":      key,
                    "old":          o_lbl,
                    "new":          n_lbl,
                    "disease":      disease,
                    "aachange":     aachange,
                    "heteroplasmy": het,
                    "homoplasmy":   hom,
                    "pubmed":       pubmed,
                    **sample_data,
                }
                for key, o_lbl, n_lbl, disease, aachange, het, hom, pubmed, sample_data
                in changed_rows
            ],
        }, f, indent=2)

    log.info("Comparacion completada")

def main():
    p = argparse.ArgumentParser(
        description="Comparador de anotaciones MitoMap (chrM)"
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("full", help="Convierte dos VCF a TSV y los compara")
    f.add_argument("old_vcf")
    f.add_argument("new_vcf")
    f.add_argument("sample_id")
    f.add_argument("--outdir", default=".")
    f.add_argument(
        "--sort-tmp-dir", default=None, metavar="DIR",
        help="Directorio temporal para el sort del sistema"
    )

    args = p.parse_args()

    if args.cmd == "full":
        out = args.outdir
        sid = args.sample_id
        old_tsv = os.path.join(out, f"{sid}.old.tsv.gz")
        new_tsv = os.path.join(out, f"{sid}.new.tsv.gz")
        vcf_to_tsv(args.old_vcf, old_tsv, sort_tmp_dir=args.sort_tmp_dir)
        vcf_to_tsv(args.new_vcf, new_tsv, sort_tmp_dir=args.sort_tmp_dir)
        compare(old_tsv, new_tsv, out, sid)


if __name__ == "__main__":
    main()