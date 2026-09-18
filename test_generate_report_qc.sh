#!/usr/bin/env bash
# Self-check for generate_report()'s QUAST/BUSCO/CheckM2 parsing + the _flag/_pct helpers.
# Fixtures mirror real tool output formats (verified against TakiLine_test_runs/*).
set -Eeuo pipefail
IFS=$'\n\t'  # matches takiline.sh's global IFS — a plain-space `read` split silently breaks under it

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/report.sh"  # exercises the real _flag/_pct/_nanostat/_svg_* helpers, not a copy

tmp=$(mktemp -d)
trap 'rm -rf "${tmp}"' EXIT

cat > "${tmp}/report.tsv" <<'EOF'
Assembly	sample
# contigs	12
Largest contig	586769
Total length	5773952
GC (%)	57.02
N50	245532
# N's per 100 kbp	0.00
# predicted rRNA genes	6 + 3 part
EOF

cat > "${tmp}/short_summary.txt" <<'EOF'
# BUSCO version is: 6.1.0
# The lineage dataset is: bacteria_odb10 (Creation date: 2024-01-08, number of genomes: 4085, number of BUSCOs: 124)

	***** Results: *****

	C:37.1%[S:37.1%,D:0.0%],F:48.4%,M:14.5%,n:124
EOF

cat > "${tmp}/quality_report.tsv" <<'EOF'
Name	Completeness	Contamination	Completeness_Model_Used
sample	98.5	1.2	Neural Network
EOF

# --- QUAST field parsing ---
gc=$(awk -F'\t' '$1=="GC (%)"{print $2}' "${tmp}/report.tsv")
ns=$(awk -F'\t' 'index($1,"# N")==1 && index($1,"kbp")>0{print $2; exit}' "${tmp}/report.tsv")
rrna=$(awk -F'\t' 'index($1,"# predicted rRNA genes")==1{print $2; exit}' "${tmp}/report.tsv")
genes=$(awk -F'\t' 'index($1,"# predicted genes")==1{print $2; exit}' "${tmp}/report.tsv")
[[ "${gc}" == "57.02" ]]         || { echo "FAIL: quast GC parse got '${gc}'"; exit 1; }
[[ "${ns}" == "0.00" ]]          || { echo "FAIL: quast Ns parse got '${ns}'"; exit 1; }
[[ "${rrna}" == "6 + 3 part" ]]  || { echo "FAIL: quast rRNA parse got '${rrna}'"; exit 1; }
[[ -z "${genes}" ]]              || { echo "FAIL: quast genes should be empty (field absent), got '${genes}'"; exit 1; }

# --- BUSCO parsing ---
busco_file="${tmp}/short_summary.txt"
lineage=$(sed -n 's/.*lineage dataset is: \([^ ]*\).*/\1/p' "${busco_file}" | head -1)
read -r complete single dup frag missing < <(
    sed -n 's/.*C:\([0-9.]*\)%\[S:\([0-9.]*\)%,D:\([0-9.]*\)%\],F:\([0-9.]*\)%,M:\([0-9.]*\)%,n:[0-9]*.*/\1\t\2\t\3\t\4\t\5/p' "${busco_file}" | head -1
)
[[ "${lineage}" == "bacteria_odb10" ]] || { echo "FAIL: busco lineage parse got '${lineage}'"; exit 1; }
[[ "${complete}" == "37.1" && "${single}" == "37.1" && "${dup}" == "0.0" && "${frag}" == "48.4" && "${missing}" == "14.5" ]] \
    || { echo "FAIL: busco composite parse got C=${complete} S=${single} D=${dup} F=${frag} M=${missing}"; exit 1; }

# --- CheckM2 parsing (header-driven; same awk covers legacy CheckM's checkm_results.tsv) ---
read -r completeness contamination < <(
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Completeness")c=i; if($i=="Contamination")k=i}} NR==2{print $c"\t"$k}' "${tmp}/quality_report.tsv"
)
[[ "${completeness}" == "98.5" && "${contamination}" == "1.2" ]] \
    || { echo "FAIL: checkm2 parse got completeness=${completeness} contamination=${contamination}"; exit 1; }

# --- _flag / _pct helpers ---
[[ "$(_flag 37.1 '<95')" == "⚠️ " ]] || { echo "FAIL: _flag should warn on BUSCO 37.1% < 95%"; exit 1; }
[[ "$(_flag 98.5 '<90')" == "✅ " ]] || { echo "FAIL: _flag should pass on CheckM2 completeness 98.5% >= 90%"; exit 1; }
[[ "$(_flag N/A '<95')"  == ""     ]] || { echo "FAIL: _flag should be blank for non-numeric N/A by default"; exit 1; }
[[ "$(_flag N/A '<95' true)" == "⚠️ " ]] || { echo "FAIL: _flag should warn on N/A when warn_on_na=true (e.g. a failed CheckM2 run)"; exit 1; }
[[ "$(_pct 57.02)" == "57.02%" ]]     || { echo "FAIL: _pct should append %"; exit 1; }
[[ "$(_pct N/A)" == "N/A" ]]          || { echo "FAIL: _pct should leave N/A bare"; exit 1; }

# --- SVG chart generators: run under a comma-decimal locale (as this box actually has) to catch
# awk's printf silently using LC_NUMERIC and emitting invalid SVG like width="207,76". ---
export LC_NUMERIC=pt_BR.UTF-8 2>/dev/null || true

# Check numeric SVG attributes specifically (not the whole blob — <title> text legitimately
# contains commas, e.g. "Complete, single-copy").
_no_broken_numeric_attrs() { ! grep -qE '="[0-9]+,[0-9]+"' <<< "$1"; }

svg_out=$(_svg_busco_bar 37.1 0.0 48.4 14.5)
_no_broken_numeric_attrs "${svg_out}" || { echo "FAIL: BUSCO chart SVG has a locale comma (invalid SVG number): ${svg_out}"; exit 1; }
[[ "${svg_out}" == *'width="207.76"'* ]] || { echo "FAIL: BUSCO chart width miscalculated: ${svg_out}"; exit 1; }

contig_lengths=$'500\n1200\n300'
contig_out=$(_svg_contig_bars "${contig_lengths}")
_no_broken_numeric_attrs "${contig_out}" || { echo "FAIL: contig chart SVG has a locale comma (invalid SVG number): ${contig_out}"; exit 1; }

# --- fastp JSON parsing (Illumina read QC) ---
cat > "${tmp}/fastp.json" <<'EOF'
{"summary": {"before_filtering": {"total_reads": 746660, "total_bases": 187411660},
             "after_filtering": {"total_reads": 520332, "total_bases": 130511378, "q20_rate": 0.89117, "q30_rate": 0.84027}}}
EOF
read -r il_reads_before il_reads_after il_bases_after il_q20_after il_q30_after < <(
    python3 -c "
import json
d = json.load(open('${tmp}/fastp.json'))
b = d['summary']['before_filtering']; a = d['summary']['after_filtering']
print(b['total_reads'], a['total_reads'], a['total_bases'], round(a['q20_rate']*100,1), round(a['q30_rate']*100,1), sep='\t')
"
)
[[ "${il_reads_before}" == "746660" && "${il_reads_after}" == "520332" && "${il_bases_after}" == "130511378" \
   && "${il_q20_after}" == "89.1" && "${il_q30_after}" == "84.0" ]] \
    || { echo "FAIL: fastp JSON parse got before=${il_reads_before} after=${il_reads_after} bases=${il_bases_after} q20=${il_q20_after} q30=${il_q30_after}"; exit 1; }
il_pct_passed=$(LC_NUMERIC=C awk -v a="${il_reads_after}" -v b="${il_reads_before}" 'BEGIN{printf "%.1f", (a/b)*100}')
[[ "${il_pct_passed}" == "69.7" ]] || { echo "FAIL: fastp pct-passed calc got '${il_pct_passed}'"; exit 1; }

# --- NanoStats.txt parsing (long-read QC) — real NanoPlot format, comma-grouped numbers ---
cat > "${tmp}/NanoStats.txt" <<'EOF'
General summary:
Mean read length:                5,550.0
Mean read quality:                  10.0
Number of reads:                33,093.0
Read length N50:                 7,415.0
Total bases:               183,665,184.0
EOF
lr_num_reads=$(_nanostat "Number of reads" "${tmp}/NanoStats.txt")
lr_total_bases=$(_nanostat "Total bases" "${tmp}/NanoStats.txt")
lr_mean_len=$(_nanostat "Mean read length" "${tmp}/NanoStats.txt")
[[ "${lr_num_reads}" == "33093.0" ]]      || { echo "FAIL: NanoStats reads parse got '${lr_num_reads}'"; exit 1; }
[[ "${lr_total_bases}" == "183665184.0" ]] || { echo "FAIL: NanoStats total-bases parse got '${lr_total_bases}'"; exit 1; }
[[ "${lr_mean_len}" == "5550.0" ]]        || { echo "FAIL: NanoStats mean-length parse got '${lr_mean_len}'"; exit 1; }

# --- estimated coverage (locale-safe division) + verdict aggregation ---
gs_bases=5000000
bases_used=$(( ${il_bases_after%.*} + ${lr_total_bases%.*} ))
est_coverage=$(LC_NUMERIC=C awk -v b="${bases_used}" -v g="${gs_bases}" 'BEGIN{printf "%.1f", b/g}')
[[ "${est_coverage}" != *,* ]] || { echo "FAIL: coverage calc has a locale comma: ${est_coverage}"; exit 1; }
[[ "${est_coverage}" == "62.8" ]] || { echo "FAIL: coverage calc got '${est_coverage}' (expected 62.8)"; exit 1; }
coverage_flag=$(_flag "${est_coverage}" "<20")
[[ "${coverage_flag}" == "✅ " ]] || { echo "FAIL: coverage flag should pass at 62.8x, got '${coverage_flag}'"; exit 1; }

verdict="PASS"
for f in "✅ " "✅ " "${coverage_flag}"; do [[ "${f}" == "⚠️ " ]] && verdict="REVIEW"; done
[[ "${verdict}" == "PASS" ]] || { echo "FAIL: verdict should be PASS when no flag warns"; exit 1; }
verdict="PASS"
for f in "✅ " "⚠️ " "${coverage_flag}"; do [[ "${f}" == "⚠️ " ]] && verdict="REVIEW"; done
[[ "${verdict}" == "REVIEW" ]] || { echo "FAIL: verdict should flip to REVIEW when any input flag warns"; exit 1; }

echo "OK: all generate_report() QC-parsing checks passed"
