import json
import math
import os
import re
import subprocess
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "outputs" / "designs"
RECEPTOR_FILE = DATA_DIR / "receptors.txt"
RECEPTOR_META_FILE = DATA_DIR / "receptor_structures.json"
RFDIFFUSION_CONFIG_FILE = DATA_DIR / "rfdiffusion_config.json"
CLOUD_CONFIG_FILE = DATA_DIR / "cloud_backend.json"
QUEUE_DIR = BASE_DIR / "outputs" / "cloud_jobs"

DEFAULT_PORT = 8000
MAX_PAPERS_LIMIT = 100
DEFAULT_PAPERS = 30
DEFAULT_BINDER_LEN = 42
MAX_BINDER_LEN = 80

EFFECT_TERMS = {
    "strong": [
        "critical",
        "essential",
        "required",
        "drives",
        "master regulator",
        "potent",
        "robust",
        "key mediator",
    ],
    "medium": [
        "promotes",
        "enhances",
        "increases",
        "activates",
        "suppresses",
        "inhibits",
        "regulates",
        "associated with",
        "linked to",
        "contributes",
    ],
}

PRO_WORDS = ["promote", "enhance", "activate", "increase", "amplify", "facilitate", "drive"]
ANTI_WORDS = ["inhibit", "suppress", "decrease", "block", "attenuate", "limit", "restrain"]

CANDIDATE_PATTERN = re.compile(
    r"\b(?:CD\d{1,3}[A-Z]?|CCR\d{1,2}|CXCR\d{1,2}|TLR\d{1,2}|IL\d{1,2}R[A-Z0-9]*|IFNGR\d?|IFNAR\d?|TNFRSF\d+[A-Z0-9]*|SIGLEC\d+|KLR[A-Z0-9]+|NCR\d|FCGR\d+[A-Z0-9]*|CSF\dR[A-Z0-9]*)\b"
)


def load_receptors():
    if not RECEPTOR_FILE.exists():
        return set()
    receptors = set()
    for line in RECEPTOR_FILE.read_text(encoding="utf-8").splitlines():
        clean = line.strip().upper()
        if clean:
            receptors.add(clean)
    return receptors


def load_receptor_meta():
    if not RECEPTOR_META_FILE.exists():
        return {}
    try:
        return json.loads(RECEPTOR_META_FILE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


KNOWN_RECEPTORS = load_receptors()
RECEPTOR_META = load_receptor_meta()

AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def safe_int(value, default, minimum=1, maximum=MAX_PAPERS_LIMIT):
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return max(minimum, min(maximum, parsed))


def safe_float(value, default, minimum=0.0, maximum=1e9):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return max(minimum, min(maximum, parsed))


def normalize_ws(text):
    return re.sub(r"\s+", " ", text or "").strip()


def split_sentences(text):
    chunks = re.split(r"(?<=[.!?])\s+", text)
    return [c.strip() for c in chunks if c.strip()]


def fetch_europe_pmc(query, page_size):
    params = {
        "query": query,
        "format": "json",
        "pageSize": page_size,
        "sort": "RELEVANCE",
    }
    encoded = urllib.parse.urlencode(params)
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{encoded}"
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "immune-receptor-explorer/1.0",
            "Accept": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    return payload.get("resultList", {}).get("result", [])


def fetch_pubmed(query, page_size):
    term = urllib.parse.urlencode({"term": query})
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=pubmed&retmode=json&retmax={page_size}&sort=relevance&{term}"
    )
    search_req = urllib.request.Request(
        search_url,
        headers={"User-Agent": "immune-receptor-explorer/1.0", "Accept": "application/json"},
    )

    with urllib.request.urlopen(search_req, timeout=30) as resp:
        search_payload = json.loads(resp.read().decode("utf-8"))

    id_list = search_payload.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return []

    id_csv = ",".join(id_list)
    fetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&retmode=xml&id={urllib.parse.quote(id_csv)}"
    )
    fetch_req = urllib.request.Request(
        fetch_url,
        headers={"User-Agent": "immune-receptor-explorer/1.0", "Accept": "application/xml"},
    )
    with urllib.request.urlopen(fetch_req, timeout=30) as resp:
        xml_text = resp.read().decode("utf-8", errors="ignore")

    root = ET.fromstring(xml_text)
    records = []
    for article in root.findall(".//PubmedArticle"):
        pmid = normalize_ws("".join(article.findtext(".//PMID", default="")))
        title = normalize_ws("".join(article.findtext(".//ArticleTitle", default="")))
        abstract_parts = [normalize_ws("".join(x.itertext())) for x in article.findall(".//Abstract/AbstractText")]
        abstract = normalize_ws(" ".join([p for p in abstract_parts if p]))
        journal = normalize_ws("".join(article.findtext(".//Journal/Title", default="")))
        year = normalize_ws("".join(article.findtext(".//PubDate/Year", default="")))
        if not year:
            medline_date = normalize_ws("".join(article.findtext(".//PubDate/MedlineDate", default="")))
            year_match = re.search(r"(19|20)\d{2}", medline_date)
            year = year_match.group(0) if year_match else ""

        records.append(
            {
                "id": pmid or "",
                "pmid": pmid or "",
                "title": title,
                "abstractText": abstract,
                "journalTitle": journal,
                "pubYear": year,
                "source": "MED",
            }
        )
    return records


def receptor_candidates(text_upper):
    found = set()
    for receptor in KNOWN_RECEPTORS:
        if re.search(rf"\b{re.escape(receptor)}\b", text_upper):
            found.add(receptor)
    for match in CANDIDATE_PATTERN.findall(text_upper):
        found.add(match)
    return found


def build_paper_link(record):
    pmid = record.get("pmid")
    if pmid:
        return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    source = record.get("source")
    ext_id = record.get("id")
    if source and ext_id:
        return f"https://europepmc.org/article/{source}/{ext_id}"
    return ""


def sentence_with_receptor(sentences, receptor):
    receptor_lower = receptor.lower()
    for s in sentences:
        if receptor_lower in s.lower():
            return s
    return ""


def score_receptors(cell, process, records):
    cell_lower = cell.lower().strip()
    process_lower = process.lower().strip()

    stats = {}

    for record in records:
        title = normalize_ws(record.get("title", ""))
        abstract = normalize_ws(record.get("abstractText", ""))
        if not (title or abstract):
            continue

        full_text = f"{title}. {abstract}".strip()
        full_upper = full_text.upper()
        sentences = split_sentences(full_text)
        candidates = receptor_candidates(full_upper)

        for receptor in candidates:
            receptor_pattern = rf"\b{re.escape(receptor)}\b"
            title_hits = len(re.findall(receptor_pattern, title.upper()))
            abstract_hits = len(re.findall(receptor_pattern, abstract.upper()))
            total_hits = title_hits + abstract_hits
            if total_hits == 0:
                continue

            receptor_stat = stats.setdefault(
                receptor,
                {
                    "receptor": receptor,
                    "title_hits": 0,
                    "abstract_hits": 0,
                    "paper_count": 0,
                    "effect_hits": 0,
                    "co_sentence_hits": 0,
                    "co_paper_hits": 0,
                    "pro_hits": 0,
                    "anti_hits": 0,
                    "papers": [],
                },
            )

            receptor_stat["title_hits"] += title_hits
            receptor_stat["abstract_hits"] += abstract_hits

            if any(p.get("paper_id") == record.get("id") for p in receptor_stat["papers"]):
                continue

            paper_has_cell = cell_lower in full_text.lower() if cell_lower else True
            paper_has_process = process_lower in full_text.lower() if process_lower else True
            if paper_has_cell and paper_has_process:
                receptor_stat["co_paper_hits"] += 1

            local_effect_hits = 0
            local_pro_hits = 0
            local_anti_hits = 0
            local_co_sentence = 0

            for s in sentences:
                s_lower = s.lower()
                if receptor.lower() not in s_lower:
                    continue

                if cell_lower and process_lower and cell_lower in s_lower and process_lower in s_lower:
                    local_co_sentence += 1
                elif cell_lower and cell_lower in s_lower:
                    local_co_sentence += 0.5
                elif process_lower and process_lower in s_lower:
                    local_co_sentence += 0.5

                for term in EFFECT_TERMS["strong"]:
                    if term in s_lower:
                        local_effect_hits += 2
                for term in EFFECT_TERMS["medium"]:
                    if term in s_lower:
                        local_effect_hits += 1

                for w in PRO_WORDS:
                    if w in s_lower:
                        local_pro_hits += 1
                for w in ANTI_WORDS:
                    if w in s_lower:
                        local_anti_hits += 1

            receptor_stat["effect_hits"] += local_effect_hits
            receptor_stat["co_sentence_hits"] += local_co_sentence
            receptor_stat["pro_hits"] += local_pro_hits
            receptor_stat["anti_hits"] += local_anti_hits
            receptor_stat["paper_count"] += 1

            snippet = sentence_with_receptor(sentences, receptor)[:400]
            receptor_stat["papers"].append(
                {
                    "paper_id": record.get("id"),
                    "pmid": record.get("pmid", ""),
                    "title": title,
                    "year": record.get("pubYear", ""),
                    "journal": record.get("journalTitle", ""),
                    "link": build_paper_link(record),
                    "snippet": snippet,
                }
            )

    results = []
    for receptor, item in stats.items():
        evidence_score = min(40, item["title_hits"] * 6 + item["abstract_hits"] * 2 + item["paper_count"] * 2)
        strength_score = min(30, item["effect_hits"] * 1.8)
        relevance_score = min(30, item["co_sentence_hits"] * 3 + item["co_paper_hits"] * 2)
        total_score = round(evidence_score + strength_score + relevance_score, 2)

        if item["pro_hits"] > item["anti_hits"] * 1.25:
            direction = "Predominantly pro-response"
        elif item["anti_hits"] > item["pro_hits"] * 1.25:
            direction = "Predominantly suppressive"
        else:
            direction = "Context-dependent / mixed"

        results.append(
            {
                "receptor": receptor,
                "total_score": total_score,
                "evidence_score": round(evidence_score, 2),
                "strength_score": round(strength_score, 2),
                "relevance_score": round(relevance_score, 2),
                "paper_count": item["paper_count"],
                "title_hits": item["title_hits"],
                "abstract_hits": item["abstract_hits"],
                "direction": direction,
                "papers": sorted(item["papers"], key=lambda p: str(p.get("year", "")), reverse=True)[:10],
            }
        )

    results.sort(key=lambda r: (r["total_score"], r["paper_count"], r["evidence_score"]), reverse=True)
    return results


def analyze(cell, process, max_papers):
    if not cell.strip() and not process.strip():
        raise ValueError("Please provide at least an immune cell or an immune process.")

    query_parts = []
    if cell.strip():
        query_parts.append(f'"{cell.strip()}"')
    if process.strip():
        query_parts.append(f'"{process.strip()}"')
    query_parts.append('("surface receptor" OR "membrane receptor" OR receptor)')
    query_parts.append('(immune OR immunology)')

    query = " AND ".join(query_parts)
    records = []
    source = ""
    errors = []

    try:
        records = fetch_europe_pmc(query, max_papers)
        source = "Europe PMC"
    except Exception as exc:
        errors.append(f"Europe PMC: {type(exc).__name__}({exc})")

    if not records:
        try:
            records = fetch_pubmed(query, max_papers)
            source = "PubMed E-utilities"
        except Exception as exc:
            errors.append(f"PubMed: {type(exc).__name__}({exc})")

    if not records and errors:
        raise RuntimeError(" | ".join(errors))

    ranked = score_receptors(cell, process, records)

    return {
        "input": {
            "cell": cell,
            "process": process,
            "max_papers": max_papers,
            "query": query,
        },
        "meta": {
            "searched_papers": len(records),
            "ranked_receptors": len(ranked),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "source": source or "Unknown",
        },
        "results": ranked[:50],
    }


def receptor_hash(receptor):
    return sum(ord(c) for c in receptor)


def generate_sequence(receptor, length):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    seed = receptor_hash(receptor)
    chars = []
    for i in range(length):
        idx = (seed + i * 7 + (i * i) * 3) % len(alphabet)
        chars.append(alphabet[idx])
    return "".join(chars)


def helix_chain(chain_id, start_resid, length, center_x, center_y, center_z):
    lines = []
    atom_id = 1
    for i in range(length):
        resid = start_resid + i
        angle = i * 1.75
        x = center_x + math.cos(angle) * 2.2
        y = center_y + math.sin(angle) * 2.2
        z = center_z + i * 1.45
        lines.append(
            f"ATOM  {atom_id:5d}  CA  ALA {chain_id}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 40.00           C"
        )
        atom_id += 1
    return lines


def build_mock_complex_pdb(receptor, binder_len):
    seed = receptor_hash(receptor)
    target_len = 120
    target_lines = []
    atom_id = 1
    for i in range(target_len):
        resid = i + 1
        theta = i * 0.32
        radius = 11.5 + ((seed + i) % 5) * 0.35
        x = math.cos(theta) * radius
        y = math.sin(theta) * radius
        z = (i - target_len / 2) * 0.75
        target_lines.append(
            f"ATOM  {atom_id:5d}  CA  GLY A{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 30.00           C"
        )
        atom_id += 1

    hotspot_center = 45 + seed % 25
    hotspot_residues = [hotspot_center - 1, hotspot_center, hotspot_center + 1, hotspot_center + 2]
    binder_center_z = (hotspot_center - target_len / 2) * 0.75

    binder_lines = []
    for i in range(binder_len):
        resid = i + 1
        angle = (seed % 9) * 0.11 + i * 1.9
        x = 14.5 + math.cos(angle) * 2.0
        y = 0.0 + math.sin(angle) * 2.0
        z = binder_center_z - 8 + i * 0.85
        binder_lines.append(
            f"ATOM  {atom_id:5d}  CA  LEU B{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
        atom_id += 1

    lines = target_lines + ["TER"] + binder_lines + ["TER", "END"]
    return "\n".join(lines) + "\n", hotspot_residues


def binding_summary(receptor, hotspot_residues):
    start = min(hotspot_residues)
    end = max(hotspot_residues)
    return {
        "target_chain": "A",
        "binder_chain": "B",
        "target_hotspot_residues": [f"A:{r}" for r in hotspot_residues],
        "binding_interface_note": f"Predicted interface around receptor residues A:{start}-A:{end}.",
        "distance_estimate_angstrom": 6.5,
        "confidence_note": "Computational prototype only. Requires structure prediction and wet-lab validation.",
    }


def receptor_icon_for_family(family):
    name = (family or "").lower()
    if "checkpoint" in name:
        return "shield"
    if "chemokine" in name:
        return "swirl"
    if "nk" in name:
        return "spark"
    return "target"


def list_target_catalog():
    catalog = []
    # The website target gallery is intentionally curated by receptor_structures.json.
    receptors = sorted(set(RECEPTOR_META.keys()))
    for receptor in receptors:
        meta = RECEPTOR_META.get(receptor, {})
        family = str(meta.get("target_family", "Unknown")).strip() or "Unknown"
        latest_dir = find_latest_design_dir(receptor)
        catalog.append(
            {
                "receptor": receptor,
                "label": meta.get("label", receptor),
                "target_family": family,
                "default_hotspot": meta.get("default_hotspot", "interface-guided"),
                "process_tags": meta.get("process_tags", []),
                "icon": receptor_icon_for_family(family),
                "has_design": bool(latest_dir),
                "latest_design_id": latest_dir.name if latest_dir else "",
            }
        )
    return catalog


def find_latest_design_dir(receptor):
    prefix = f"{receptor.strip().upper()}_"
    candidates = [p for p in OUTPUT_DIR.glob(f"{prefix}*") if p.is_dir()]
    if not candidates:
        return None
    candidates.sort(key=lambda p: p.name, reverse=True)
    return candidates[0]


def load_design_payload_from_dir(design_dir, receptor):
    report_file = design_dir / "design_report.json"
    report = {}
    if report_file.exists():
        try:
            report = json.loads(report_file.read_text(encoding="utf-8"))
            if not isinstance(report, dict):
                report = {}
        except json.JSONDecodeError:
            report = {}

    complex_file = design_dir / "complex_model.pdb"
    binder_file = design_dir / "binder_chain_B.pdb"
    if not complex_file.exists():
        pdb_candidates = list_output_pdb_files(design_dir)
        auto_complex, auto_binder = choose_complex_and_binder_pdb(pdb_candidates)
        if auto_complex:
            complex_file = auto_complex
        if auto_binder:
            binder_file = auto_binder

    if not complex_file.exists():
        raise FileNotFoundError("No complex model PDB was found for this design.")

    complex_text = complex_file.read_text(encoding="utf-8", errors="ignore")
    binder_text = ""
    if binder_file.exists():
        binder_text = binder_file.read_text(encoding="utf-8", errors="ignore")
    else:
        tmp_binder = design_dir / "binder_chain_B.pdb"
        if write_chain_b_from_complex(complex_text, tmp_binder):
            binder_file = tmp_binder
            binder_text = tmp_binder.read_text(encoding="utf-8", errors="ignore")

    design_id = str(report.get("design_id", design_dir.name))
    receptor_out = str(report.get("receptor", receptor))
    binder_seq = str(report.get("binder_sequence", "")).strip()
    if not binder_seq:
        binder_seq = parse_pdb_chain_sequence(complex_text, "B") or "N/A"

    files_payload = {
        "complex_pdb": f"/designs/{design_id}/{complex_file.name}",
        "binder_pdb": f"/designs/{design_id}/{binder_file.name}" if binder_file.exists() else "",
        "report_json": f"/designs/{design_id}/design_report.json",
    }

    payload = {
        "design_id": design_id,
        "receptor": receptor_out,
        "cell": report.get("cell", ""),
        "process": report.get("process", ""),
        "mode": report.get("mode", "prototype"),
        "engine": report.get("engine", "loaded-from-disk"),
        "created_at": report.get("created_at", ""),
        "binder_length": report.get("binder_length", len(binder_seq) if binder_seq != "N/A" else 0),
        "binder_sequence": binder_seq,
        "receptor_annotation": report.get("receptor_annotation", RECEPTOR_META.get(receptor_out, {})),
        "requested_hotspot": report.get("requested_hotspot", ""),
        "binding": report.get(
            "binding",
            {
                "target_chain": "A",
                "binder_chain": "B",
                "target_hotspot_residues": [],
                "binding_interface_note": "Loaded existing design. Interface note not available.",
                "distance_estimate_angstrom": "n/a",
                "confidence_note": "Loaded from local design artifacts.",
            },
        ),
        "files": files_payload,
        "warnings": report.get("warnings", []),
        "complex_pdb_text": complex_text,
        "binder_pdb_text": binder_text,
    }
    return payload


def get_latest_design_for_receptor(receptor):
    receptor = receptor.strip().upper()
    if not receptor:
        raise ValueError("receptor is required")
    design_dir = find_latest_design_dir(receptor)
    if not design_dir:
        return None
    return load_design_payload_from_dir(design_dir, receptor)


def load_rfdiffusion_config():
    if not RFDIFFUSION_CONFIG_FILE.exists():
        return {}
    try:
        payload = json.loads(RFDIFFUSION_CONFIG_FILE.read_text(encoding="utf-8"))
        return payload if isinstance(payload, dict) else {}
    except json.JSONDecodeError:
        return {}


def load_cloud_config():
    if not CLOUD_CONFIG_FILE.exists():
        return {}
    try:
        payload = json.loads(CLOUD_CONFIG_FILE.read_text(encoding="utf-8"))
        return payload if isinstance(payload, dict) else {}
    except json.JSONDecodeError:
        return {}


def cloud_status():
    cfg = load_cloud_config()
    submit_url = str(cfg.get("submit_url", "")).strip()
    status_url = str(cfg.get("status_url_template", "")).strip()
    result_url = str(cfg.get("result_url_template", "")).strip()
    configured = bool(submit_url and status_url and result_url)
    return {
        "configured": configured,
        "config_file": str(CLOUD_CONFIG_FILE),
        "submit_url_set": bool(submit_url),
        "status_url_set": bool(status_url),
        "result_url_set": bool(result_url),
        "note": "Cloud mode requires submit_url, status_url_template and result_url_template.",
    }


def rfdiffusion_status():
    cfg = load_rfdiffusion_config()
    env_cmd = os.environ.get("RFDIFFUSION_CMD", "").strip()
    configured = bool(cfg.get("command") or env_cmd)
    return {
        "configured": configured,
        "config_file": str(RFDIFFUSION_CONFIG_FILE),
        "uses_env_command": bool(env_cmd),
        "uses_config_command": bool(cfg.get("command")),
        "engine_ready_note": "Configured means command is provided. It does not guarantee runtime success.",
    }


def http_json_request(url, method="GET", payload=None, headers=None, timeout=30):
    hdrs = {"Accept": "application/json", "User-Agent": "immune-receptor-explorer/1.0"}
    if headers:
        hdrs.update(headers)
    data = None
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        hdrs["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8", errors="ignore")
    parsed = json.loads(body) if body.strip() else {}
    return parsed if isinstance(parsed, dict) else {"raw": parsed}


def write_cloud_job_audit(design_id, payload):
    QUEUE_DIR.mkdir(parents=True, exist_ok=True)
    out = QUEUE_DIR / f"{design_id}.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return out


def submit_cloud_design_job(receptor, cell, process, binder_len, hotspot_text, design_id):
    cfg = load_cloud_config()
    submit_url = str(cfg.get("submit_url", "")).strip()
    if not submit_url:
        raise RuntimeError("Cloud backend is not configured. Please set data/cloud_backend.json.")

    api_key_env = str(cfg.get("api_key_env", "CLOUD_RF_API_KEY")).strip() or "CLOUD_RF_API_KEY"
    api_key = os.environ.get(api_key_env, "").strip()
    timeout_sec = safe_int(cfg.get("timeout_sec", 60), 60, minimum=10, maximum=600)
    headers = {}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    payload = {
        "design_id": design_id,
        "receptor": receptor,
        "cell": cell,
        "process": process,
        "binder_length": binder_len,
        "hotspot_hint": hotspot_text or "",
        "callback_expected": False,
    }
    response = http_json_request(submit_url, method="POST", payload=payload, headers=headers, timeout=timeout_sec)
    job_id = str(response.get("job_id", design_id)).strip() or design_id
    status = str(response.get("status", "queued")).strip() or "queued"
    audit = {
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "submit_url": submit_url,
        "request": payload,
        "response": response,
    }
    write_cloud_job_audit(design_id, audit)
    return {
        "status": status,
        "job_id": job_id,
        "design_id": design_id,
        "engine": "rfdiffusion-cloud",
        "message": "Cloud job submitted. Poll status until finished.",
        "cloud_poll_url": f"/api/cloud_job?job_id={urllib.parse.quote(job_id)}",
    }


def poll_cloud_job(job_id):
    cfg = load_cloud_config()
    status_tpl = str(cfg.get("status_url_template", "")).strip()
    result_tpl = str(cfg.get("result_url_template", "")).strip()
    if not status_tpl or not result_tpl:
        raise RuntimeError("Cloud backend status/result URL template is not configured.")

    api_key_env = str(cfg.get("api_key_env", "CLOUD_RF_API_KEY")).strip() or "CLOUD_RF_API_KEY"
    api_key = os.environ.get(api_key_env, "").strip()
    timeout_sec = safe_int(cfg.get("timeout_sec", 60), 60, minimum=10, maximum=600)
    headers = {}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    status_url = status_tpl.format(job_id=job_id)
    status_payload = http_json_request(status_url, method="GET", payload=None, headers=headers, timeout=timeout_sec)
    status = str(status_payload.get("status", "")).lower().strip() or "unknown"
    if status not in {"done", "finished", "completed", "success"}:
        return {
            "job_id": job_id,
            "status": status,
            "engine": "rfdiffusion-cloud",
            "message": status_payload.get("message", "Cloud job still running."),
            "status_payload": status_payload,
        }

    result_url = result_tpl.format(job_id=job_id)
    result_payload = http_json_request(result_url, method="GET", payload=None, headers=headers, timeout=timeout_sec)
    if "design_id" not in result_payload:
        result_payload["design_id"] = job_id
    result_payload["engine"] = result_payload.get("engine", "rfdiffusion-cloud")
    result_payload["mode"] = result_payload.get("mode", "rfdiffusion-cloud")
    result_payload["status"] = "completed"
    return result_payload


def parse_pdb_chain_sequence(pdb_text, chain_id):
    seen = set()
    seq = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 26:
            continue
        chain = line[21:22]
        if chain != chain_id:
            continue
        resname = line[17:20].strip().upper()
        resid = line[22:26].strip()
        key = (chain, resid)
        if key in seen:
            continue
        seen.add(key)
        seq.append(AA3_TO_AA1.get(resname, "X"))
    return "".join(seq)


def list_output_pdb_files(out_dir):
    if not out_dir.exists():
        return []
    return sorted(out_dir.glob("*.pdb"), key=lambda p: p.stat().st_size, reverse=True)


def choose_complex_and_binder_pdb(pdb_files):
    if not pdb_files:
        return None, None
    complex_file = pdb_files[0]
    binder_file = None
    for f in pdb_files:
        name = f.name.lower()
        if "binder" in name or "_b" in name or "chainb" in name:
            binder_file = f
            break
    return complex_file, binder_file


def write_chain_b_from_complex(complex_text, out_path):
    chain_lines = []
    for line in complex_text.splitlines():
        if line.startswith("ATOM") and len(line) > 21 and line[21:22] == "B":
            chain_lines.append(line)
    if not chain_lines:
        return False
    chain_lines.extend(["TER", "END"])
    out_path.write_text("\n".join(chain_lines) + "\n", encoding="utf-8")
    return True


def render_command_template(template, context):
    rendered = []
    for part in template:
        rendered.append(str(part).format(**context))
    return rendered


def maybe_run_rfdiffusion(receptor, binder_len, hotspot_text, cell, process, design_id, out_dir):
    cfg = load_rfdiffusion_config()
    env_cmd = os.environ.get("RFDIFFUSION_CMD", "").strip()

    command = cfg.get("command")
    if not command and env_cmd:
        command = [env_cmd, "{receptor}", "{binder_len}", "{hotspot_hint}", "{output_prefix}"]

    if not command:
        return None
    if not isinstance(command, list):
        raise RuntimeError("rfdiffusion_config.command must be a JSON array of command tokens.")

    timeout_sec = safe_int(cfg.get("timeout_sec", 1800), 1800, minimum=60, maximum=7200)
    workdir = str(Path(cfg.get("workdir", ".")).resolve())
    output_prefix = str((out_dir / design_id).resolve())

    context = {
        "receptor": receptor,
        "binder_len": binder_len,
        "hotspot_hint": hotspot_text or "",
        "cell": cell or "",
        "process": process or "",
        "design_id": design_id,
        "out_dir": str(out_dir.resolve()),
        "output_prefix": output_prefix,
    }
    cmd = render_command_template(command, context)

    completed = subprocess.run(
        cmd,
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=timeout_sec,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"RFdiffusion command failed (code={completed.returncode}). stderr={completed.stderr[-5000:]}"
        )

    stdout_text = (completed.stdout or "").strip()
    if stdout_text.startswith("{") and stdout_text.endswith("}"):
        try:
            data = json.loads(stdout_text)
            if isinstance(data, dict):
                data["engine"] = "rfdiffusion"
                data["mode"] = "rfdiffusion"
                return data
        except json.JSONDecodeError:
            pass

    pdb_files = list_output_pdb_files(out_dir)
    complex_file, binder_file = choose_complex_and_binder_pdb(pdb_files)
    if not complex_file:
        raise RuntimeError("RFdiffusion finished but no PDB files were found in output directory.")

    complex_text = complex_file.read_text(encoding="utf-8", errors="ignore")
    if not binder_file:
        binder_file = out_dir / "binder_chain_B.pdb"
        write_chain_b_from_complex(complex_text, binder_file)

    binder_text = ""
    if binder_file.exists():
        binder_text = binder_file.read_text(encoding="utf-8", errors="ignore")

    binder_seq = parse_pdb_chain_sequence(complex_text, "B")
    if not binder_seq and binder_text:
        binder_seq = parse_pdb_chain_sequence(binder_text, "B")

    hotspot_res = [f"A:{r}" for r in range(60, 64)]
    return {
        "design_id": design_id,
        "receptor": receptor,
        "cell": cell,
        "process": process,
        "mode": "rfdiffusion",
        "engine": "rfdiffusion",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "binder_length": len(binder_seq) if binder_seq else binder_len,
        "binder_sequence": binder_seq or "N/A",
        "receptor_annotation": RECEPTOR_META.get(receptor, {}),
        "requested_hotspot": hotspot_text or "interface-guided",
        "binding": {
            "target_chain": "A",
            "binder_chain": "B",
            "target_hotspot_residues": hotspot_res,
            "binding_interface_note": "Interface residues should be refined using docking/AF2 confidence metrics.",
            "distance_estimate_angstrom": 6.0,
            "confidence_note": "Generated by RFdiffusion pipeline. Experimental validation required.",
        },
        "files": {
            "complex_pdb": f"/designs/{design_id}/{complex_file.name}",
            "binder_pdb": f"/designs/{design_id}/{binder_file.name}" if binder_file else "",
            "report_json": f"/designs/{design_id}/design_report.json",
        },
        "warnings": [
            "RFdiffusion output is computational and must be validated experimentally.",
        ],
        "complex_pdb_text": complex_text,
        "binder_pdb_text": binder_text,
        "rf_run": {
            "command": cmd,
            "workdir": workdir,
        },
    }


def design_mini_binder(receptor, cell, process, binder_len, hotspot_text="", execution_mode="mock"):
    receptor = receptor.strip().upper()
    if not receptor:
        raise ValueError("receptor is required")

    binder_len = safe_int(binder_len, DEFAULT_BINDER_LEN, minimum=20, maximum=MAX_BINDER_LEN)
    now = datetime.now(timezone.utc)
    design_id = f"{receptor}_{now.strftime('%Y%m%dT%H%M%SZ')}"
    out_dir = OUTPUT_DIR / design_id
    out_dir.mkdir(parents=True, exist_ok=True)

    mode = (execution_mode or "mock").strip().lower()
    if mode == "cloud":
        return submit_cloud_design_job(receptor, cell, process, binder_len, hotspot_text, design_id)
    if mode == "mock":
        rfd_result = None
    else:
        rfd_result = maybe_run_rfdiffusion(receptor, binder_len, hotspot_text, cell, process, design_id, out_dir)
    if rfd_result:
        report_file = out_dir / "design_report.json"
        report_file.write_text(json.dumps(rfd_result, ensure_ascii=False, indent=2), encoding="utf-8")
        return rfd_result

    receptor_info = RECEPTOR_META.get(receptor, {})
    sequence = generate_sequence(receptor, binder_len)
    complex_pdb, hotspot_residues = build_mock_complex_pdb(receptor, binder_len)

    complex_file = out_dir / "complex_model.pdb"
    binder_file = out_dir / "binder_chain_B.pdb"
    report_file = out_dir / "design_report.json"

    complex_file.write_text(complex_pdb, encoding="utf-8")

    binder_only = []
    for line in complex_pdb.splitlines():
        if line.startswith("ATOM") and line[21:22] == "B":
            binder_only.append(line)
    binder_only.append("TER")
    binder_only.append("END")
    binder_file.write_text("\n".join(binder_only) + "\n", encoding="utf-8")

    report = {
        "design_id": design_id,
        "receptor": receptor,
        "cell": cell,
        "process": process,
        "mode": "prototype",
        "engine": "mock-rfdiffusion-framework",
        "created_at": now.isoformat(),
        "binder_length": binder_len,
        "binder_sequence": sequence,
        "receptor_annotation": receptor_info,
        "requested_hotspot": hotspot_text or receptor_info.get("default_hotspot", "interface-guided"),
        "binding": binding_summary(receptor, hotspot_residues),
        "files": {
            "complex_pdb": f"/designs/{design_id}/complex_model.pdb",
            "binder_pdb": f"/designs/{design_id}/binder_chain_B.pdb",
            "report_json": f"/designs/{design_id}/design_report.json",
        },
        "warnings": [
            "This output is a computational prototype for ideation.",
            "It is not a validated therapeutic candidate.",
            "Run full RFdiffusion + AF2/ESMFold + docking + wet-lab assays before conclusions.",
        ],
    }
    report_file.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    report["complex_pdb_text"] = complex_pdb
    report["binder_pdb_text"] = binder_file.read_text(encoding="utf-8")
    return report


class AppHandler(BaseHTTPRequestHandler):
    def send_json(self, payload, code=200):
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def parse_json_body(self):
        length = int(self.headers.get("Content-Length", "0"))
        if length <= 0:
            return {}
        raw = self.rfile.read(length)
        return json.loads(raw.decode("utf-8"))

    def send_file(self, path):
        if not path.exists() or not path.is_file():
            self.send_error(404, "File not found")
            return

        suffix = path.suffix.lower()
        mime_types = {
            ".html": "text/html; charset=utf-8",
            ".css": "text/css; charset=utf-8",
            ".js": "application/javascript; charset=utf-8",
            ".json": "application/json; charset=utf-8",
            ".txt": "text/plain; charset=utf-8",
            ".pdb": "chemical/x-pdb; charset=utf-8",
        }
        content_type = mime_types.get(suffix, "application/octet-stream")
        content = path.read_bytes()

        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(content)))
        self.end_headers()
        self.wfile.write(content)

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        route = parsed.path

        if route == "/api/health":
            self.send_json({"status": "ok", "service": "immune-receptor-explorer"})
            return
        if route == "/api/rfdiffusion_status":
            self.send_json(rfdiffusion_status())
            return
        if route == "/api/cloud_status":
            self.send_json(cloud_status())
            return
        if route == "/api/cloud_job":
            params = urllib.parse.parse_qs(parsed.query)
            job_id = params.get("job_id", [""])[0].strip()
            if not job_id:
                self.send_json({"error": "job_id is required"}, code=400)
                return
            try:
                payload = poll_cloud_job(job_id)
                self.send_json(payload)
            except Exception as exc:
                self.send_json({"error": "Cloud polling failed.", "detail": repr(exc)}, code=502)
            return

        if route == "/api/targets":
            self.send_json({"targets": list_target_catalog()})
            return

        if route == "/api/target_design":
            params = urllib.parse.parse_qs(parsed.query)
            receptor = params.get("receptor", [""])[0].strip().upper()
            if not receptor:
                self.send_json({"error": "receptor is required"}, code=400)
                return
            try:
                result = get_latest_design_for_receptor(receptor)
            except Exception as exc:
                self.send_json({"error": "Failed to load latest design.", "detail": repr(exc)}, code=500)
                return
            if not result:
                self.send_json({"status": "not_found", "receptor": receptor}, code=404)
                return
            self.send_json(result)
            return

        if route == "/api/analyze":
            params = urllib.parse.parse_qs(parsed.query)
            cell = params.get("cell", [""])[0]
            process = params.get("process", [""])[0]
            max_papers = safe_int(params.get("max_papers", [DEFAULT_PAPERS])[0], DEFAULT_PAPERS)

            try:
                payload = analyze(cell, process, max_papers)
                self.send_json(payload)
            except ValueError as exc:
                self.send_json({"error": str(exc)}, code=400)
            except Exception as exc:
                self.send_json(
                    {
                        "error": "Failed to query literature databases.",
                        "detail": repr(exc),
                    },
                    code=502,
                )
            return

        if route.startswith("/designs/"):
            target = (OUTPUT_DIR.parent / route.lstrip("/")).resolve()
            if OUTPUT_DIR.parent.resolve() not in target.parents and target != OUTPUT_DIR.parent.resolve():
                self.send_error(403, "Forbidden")
                return
            self.send_file(target)
            return

        if route in ["/", ""]:
            self.send_file(STATIC_DIR / "index.html")
            return

        static_path = (STATIC_DIR / route.lstrip("/")).resolve()
        if STATIC_DIR.resolve() not in static_path.parents and static_path != STATIC_DIR.resolve():
            self.send_error(403, "Forbidden")
            return

        self.send_file(static_path)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path != "/api/design_mini_binder":
            self.send_error(404, "Not found")
            return

        try:
            payload = self.parse_json_body()
            receptor = str(payload.get("receptor", ""))
            cell = str(payload.get("cell", ""))
            process = str(payload.get("process", ""))
            binder_len = payload.get("binder_length", DEFAULT_BINDER_LEN)
            hotspot_text = str(payload.get("hotspot_hint", ""))
            execution_mode = str(payload.get("execution_mode", "mock"))

            result = design_mini_binder(receptor, cell, process, binder_len, hotspot_text, execution_mode)
            self.send_json(result)
        except ValueError as exc:
            self.send_json({"error": str(exc)}, code=400)
        except Exception as exc:
            self.send_json({"error": "Mini-binder design failed.", "detail": repr(exc)}, code=500)


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    QUEUE_DIR.mkdir(parents=True, exist_ok=True)
    host = os.environ.get("HOST", "0.0.0.0").strip() or "0.0.0.0"
    raw_port = os.environ.get("PORT", str(DEFAULT_PORT)).strip()
    try:
        port = int(raw_port)
    except ValueError:
        port = DEFAULT_PORT
    server = ThreadingHTTPServer((host, port), AppHandler)
    print(f"Server running on http://{host}:{port}")
    server.serve_forever()
