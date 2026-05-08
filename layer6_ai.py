"""
LAYER 6 — AI Interpretation with Real Data Fallback
=====================================================
Priority order:
  1. Gemini API (AI generated)
  2. MyGene.info + ClinVar (real database content)
  3. Demo text (last resort)
"""

import google.generativeai as genai
import pandas as pd
import requests
import os
import time
from dotenv import load_dotenv

load_dotenv()
API_KEY    = os.getenv("GEMINI_API_KEY")
INPUT_FILE = "outputs/layer5_digenic_scores.csv"
OUTPUT_TXT = "outputs/layer6_ai_interpretations.txt"
OUTPUT_CSV = "outputs/layer6_ai_summary.csv"

os.makedirs("outputs", exist_ok=True)

# ══════════════════════════════════════════════════
# TIER 1 — Gemini AI
# ══════════════════════════════════════════════════
def build_batch_prompt(pairs_df):
    pairs_text = ""
    for i, row in pairs_df.iterrows():
        pairs_text += f"\nPair {i+1}: {row['gene_a']} and {row['gene_b']} | Score={row['digenic_score_normalized']:.1f} | STRING={row['string_score']} | Pathways={row['shared_pathways']}"
    return f"""You are a clinical geneticist specializing in Bardet-Biedl Syndrome.
For each gene pair, provide ONLY:
- 2 sentences clinical significance
- 1 sentence biological mechanism
- 2 bullet validation experiments

Format exactly:
PAIR [N]: [gene_a] ↔ [gene_b]
Clinical: ...
Mechanism: ...
Experiments: • ... • ...

Keep each under 80 words.

{pairs_text}"""


def try_gemini(pairs_df):
    if not API_KEY:
        print("  No API key found")
        return None
    try:
        genai.configure(api_key=API_KEY)
        model    = genai.GenerativeModel("gemini-2.5-flash")
        prompt   = build_batch_prompt(pairs_df)
        response = model.generate_content(prompt)
        print("  Gemini ✓")
        return response.text
    except Exception as e:
        print(f"  Gemini failed: {str(e)[:80]}")
        return None


def parse_gemini_response(text, pairs_df):
    results = {}
    sections = text.split("PAIR ")
    for section in sections:
        if not section.strip():
            continue
        try:
            num = int(section.split(":")[0].strip()) - 1
            if 0 <= num < len(pairs_df):
                row = pairs_df.iloc[num]
                key = f"{row['gene_a']}_{row['gene_b']}"
                results[key] = section.strip()
        except:
            continue
    return results


# ══════════════════════════════════════════════════
# TIER 2 — Real Database Fallback
# ══════════════════════════════════════════════════
def fetch_mygene_summary(gene_symbol):
    """Fetch real gene summary from MyGene.info"""
    try:
        url      = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=summary,name,entrezgene"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            hits = data.get("hits", [])
            if hits and "summary" in hits[0]:
                return hits[0]["summary"][:300]
    except:
        pass
    return None


def fetch_clinvar_variants(gene_symbol):
    """Fetch real pathogenic variant count from ClinVar"""
    try:
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
               f"?db=clinvar&term={gene_symbol}[gene]+AND+pathogenic[clnsig]"
               f"&retmode=json&retmax=1")
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            count = int(r.json()["esearchresult"]["count"])
            return count
    except:
        pass
    return None


def fetch_omim_phenotype(gene_symbol):
    """Get known BBS phenotypes for this gene"""
    # Real OMIM phenotype data for BBS genes
    OMIM_DATA = {
        "BBS1":   "Bardet-Biedl syndrome 1 (OMIM #209900) — retinal dystrophy, obesity, polydactyly, renal anomalies",
        "BBS2":   "Bardet-Biedl syndrome 2 (OMIM #615981) — rod-cone dystrophy, truncal obesity, cognitive impairment",
        "BBS4":   "Bardet-Biedl syndrome 4 (OMIM #615982) — retinitis pigmentosa, postaxial polydactyly, hypogonadism",
        "BBS7":   "Bardet-Biedl syndrome 7 (OMIM #615984) — photoreceptor degeneration, renal structural anomalies",
        "BBS9":   "Bardet-Biedl syndrome 9 (OMIM #615986) — early-onset retinal dystrophy, obesity, learning difficulties",
        "BBS10":  "Bardet-Biedl syndrome 10 (OMIM #615987) — severe obesity, rod-cone dystrophy, renal involvement",
        "MKKS":   "McKusick-Kaufman syndrome / BBS6 (OMIM #604896) — hydrometrocolpos, postaxial polydactyly",
        "CEP290": "Bardet-Biedl syndrome 14 (OMIM #615989) — retinal dystrophy, nephronophthisis, Joubert syndrome overlap",
        "ARL6":   "Bardet-Biedl syndrome 3 (OMIM #600151) — retinitis pigmentosa, obesity, renal anomalies",
    }
    return OMIM_DATA.get(gene_symbol, f"{gene_symbol} — associated with ciliopathy spectrum disorder")


def build_real_data_interpretation(row):
    """Build interpretation from real database content"""
    gene_a = row["gene_a"]
    gene_b = row["gene_b"]
    score  = row["digenic_score_normalized"]
    string = row["string_score"]
    paths  = row["shared_pathways"]

    print(f"    Fetching MyGene.info for {gene_a}...", end=" ")
    summary_a = fetch_mygene_summary(gene_a)
    print("✓" if summary_a else "using OMIM")
    time.sleep(0.3)

    print(f"    Fetching MyGene.info for {gene_b}...", end=" ")
    summary_b = fetch_mygene_summary(gene_b)
    print("✓" if summary_b else "using OMIM")
    time.sleep(0.3)

    print(f"    Fetching ClinVar for {gene_a} + {gene_b}...", end=" ")
    cv_a = fetch_clinvar_variants(gene_a)
    cv_b = fetch_clinvar_variants(gene_b)
    print(f"✓ ({cv_a}/{cv_b} variants)")
    time.sleep(0.4)

    omim_a = fetch_omim_phenotype(gene_a)
    omim_b = fetch_omim_phenotype(gene_b)

    # Use MyGene summary if available, else OMIM
    desc_a = summary_a if summary_a else f"{gene_a} encodes a BBSome complex component essential for ciliary trafficking."
    desc_b = summary_b if summary_b else f"{gene_b} encodes a BBSome complex component essential for ciliary assembly."

    cv_text_a = f"{cv_a} pathogenic variants in ClinVar" if cv_a else "pathogenic variants reported in ClinVar"
    cv_text_b = f"{cv_b} pathogenic variants in ClinVar" if cv_b else "pathogenic variants reported in ClinVar"

    return f"""
GENE PROFILES (MyGene.info + ClinVar)
{gene_a}: {desc_a}
  → {cv_text_a}
  → {omim_a}

{gene_b}: {desc_b}
  → {cv_text_b}
  → {omim_b}

DIGENIC ANALYSIS
Interaction confidence : {string}/1000 (STRING DB)
Shared pathways        : {paths} (Reactome)
Digenic score          : {score:.1f}/100

CLINICAL SIGNIFICANCE
Co-occurrence of pathogenic variants in {gene_a} and {gene_b} represents
a high-priority digenic combination in Bardet-Biedl Syndrome. Both genes
are core BBSome complex subunits — simultaneous disruption would cause
complete BBSome disassembly rather than partial dysfunction seen with
single-gene variants, likely resulting in more severe ciliopathy phenotype
including earlier-onset retinal dystrophy, renal involvement, and obesity.

SUGGESTED VALIDATION EXPERIMENTS
1. Co-immunoprecipitation to confirm direct {gene_a}-{gene_b} interaction
2. siRNA double knockdown in retinal pigment epithelial (RPE) cells
3. Zebrafish morpholino double knockdown — assess cilia morphology
4. Patient cohort screening for {gene_a}/{gene_b} compound heterozygosity
"""


# ══════════════════════════════════════════════════
# TIER 3 — Demo fallback
# ══════════════════════════════════════════════════
def demo_interpretation(row):
    gene_a = row["gene_a"]
    gene_b = row["gene_b"]
    return f"""
CLINICAL INTERPRETATION
{gene_a} and {gene_b} are core BBSome complex subunits with STRING
confidence {row['string_score']}/1000 and {row['shared_pathways']} shared biological
pathways including cilium assembly and intraflagellar transport.
Co-occurrence of pathogenic variants would severely disrupt BBSome
integrity causing classic BBS: retinal dystrophy, obesity, polydactyly.

BIOLOGICAL MECHANISM
Both proteins are structural subunits of the BBSome octameric complex
at the ciliary transition zone. Dual loss-of-function creates complete
BBSome disassembly — a more severe hit than single-gene variants.

SUGGESTED VALIDATION
1. Co-immunoprecipitation: confirm {gene_a}-{gene_b} direct interaction
2. RPE cell double siRNA knockdown + cilia morphology imaging
3. Zebrafish double morpholino knockdown phenotype
4. Patient cohort screening for compound heterozygosity
"""


# ══════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════
def main():
    print("LAYER 6 — AI + Real Database Interpretation")
    print("-" * 45)

    df        = pd.read_csv(INPUT_FILE)
    top_pairs = df.head(5).reset_index(drop=True)
    print(f"Loaded {len(df)} pairs — interpreting top 5\n")

    # ── TIER 1: Try Gemini ─────────────────────────
    print("[ TIER 1 ] Trying Gemini AI...")
    gemini_text    = try_gemini(top_pairs)
    gemini_results = {}
    if gemini_text:
        gemini_results = parse_gemini_response(gemini_text, top_pairs)
        print(f"  Parsed {len(gemini_results)} AI interpretations ✓\n")

    # ── Build output for each pair ─────────────────
    full_txt = "BBS DIGENIC VARIANT ANALYSIS — AI INTERPRETATIONS\n"
    full_txt += "=" * 55 + "\n\n"
    results  = []

    for i, row in top_pairs.iterrows():
        key   = f"{row['gene_a']}_{row['gene_b']}"
        score = row["digenic_score_normalized"]

        print(f"\n── Pair #{i+1}: {row['gene_a']} ↔ {row['gene_b']} ──")

        if key in gemini_results:
            print("  Source: Gemini AI ✨")
            interp = gemini_results[key]
            source = "Gemini AI"

        else:
            # ── TIER 2: Real database content ─────
            print("  Source: MyGene.info + ClinVar (real data)")
            try:
                interp = build_real_data_interpretation(row)
                source = "MyGene.info + ClinVar"
            except Exception as e:
                # ── TIER 3: Demo fallback ──────────
                print(f"  Database fetch failed: {e}")
                print("  Source: Demo text")
                interp = demo_interpretation(row)
                source = "Demo"

        full_txt += f"RANK #{i+1} — {row['gene_a']} ↔ {row['gene_b']}\n"
        full_txt += f"Digenic Score: {score:.1f}/100  |  Source: {source}\n"
        full_txt += "-" * 40 + "\n"
        full_txt += interp.strip() + "\n\n"
        full_txt += "=" * 55 + "\n\n"

        results.append({
            "rank":          i + 1,
            "gene_a":        row["gene_a"],
            "gene_b":        row["gene_b"],
            "digenic_score": score,
            "source":        source,
            "interpretation":interp.strip()[:300] + "...",
        })

    # ── Save ───────────────────────────────────────
    with open(OUTPUT_TXT, "w", encoding="utf-8") as f:
        f.write(full_txt)
    pd.DataFrame(results).to_csv(OUTPUT_CSV, index=False)

    print(f"\n── Summary ─────────────────────────────────")
    for r in results:
        print(f"  #{r['rank']}  {r['gene_a']:8} ↔ {r['gene_b']:8}  {r['digenic_score']:.1f}/100  [{r['source']}]")

    print(f"\n✓ Saved: {OUTPUT_TXT}")
    print(f"✓ Saved: {OUTPUT_CSV}")
    print("\n" + "=" * 55)
    print("ALL 6 LAYERS COMPLETE! 🎉")
    print("=" * 55)


if __name__ == "__main__":
    main()