_flag() {
    local v="$1" cond="$2" warn_on_na="${3:-false}"
    if [[ ! "${v}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        ${warn_on_na} && echo -n "⚠️ "
        return 0
    fi
    awk "BEGIN{exit !(${v}${cond})}" && echo -n "⚠️ " || echo -n "✅ "
}

_pct() { [[ "$1" == "N/A" ]] && echo -n "N/A" || echo -n "${1}%"; }

_nanostat() { awk -F: -v k="$1" '$1==k{v=$2; gsub(/[ ,]/,"",v); print v}' "$2"; }

_svg_busco_bar() {
    local single="$1" dup="$2" frag="$3" missing="$4" v
    for v in "${single}" "${dup}" "${frag}" "${missing}"; do
        [[ "${v}" =~ ^[0-9]+(\.[0-9]+)?$ ]] || return 0
    done
    LC_NUMERIC=C awk -v s="${single}" -v d="${dup}" -v f="${frag}" -v m="${missing}" '
        BEGIN {
            W = 560; H = 36; xs = 0
            printf "<rect x=\"0\" y=\"0\" width=\"%.2f\" height=\"%d\" fill=\"#e99e44\"><title>Complete, single-copy: %s%%</title></rect>\n", (s/100)*W, H, s
            xs += (s/100)*W
            printf "<rect x=\"%.2f\" y=\"0\" width=\"%.2f\" height=\"%d\" fill=\"#3492c7\"><title>Complete, duplicated: %s%%</title></rect>\n", xs, (d/100)*W, H, d
            xs += (d/100)*W
            printf "<rect x=\"%.2f\" y=\"0\" width=\"%.2f\" height=\"%d\" fill=\"#f0c419\"><title>Fragmented: %s%%</title></rect>\n", xs, (f/100)*W, H, f
            xs += (f/100)*W
            printf "<rect x=\"%.2f\" y=\"0\" width=\"%.2f\" height=\"%d\" fill=\"#c0392b\"><title>Missing: %s%%</title></rect>\n", xs, (m/100)*W, H, m
        }'
}

_svg_contig_chart() {
    local lengths="$1"
    [[ -n "${lengths}" ]] || return 0
    sort -rn <<< "${lengths}" | LC_NUMERIC=C awk '
        BEGIN { W=700; H=260; PAD_L=75; PAD_R=15; PAD_T=10; PAD_B=26; PW=W-PAD_L-PAD_R; PH=H-PAD_T-PAD_B }
        { len[NR]=$1; total+=$1; n=NR }
        END {
            if (n==0 || total==0) exit
            maxlen=len[1]; cum=0; n50_x=-1
            path=sprintf("M %.2f,%.2f", PAD_L, PAD_T+PH)
            for (i=1; i<=n; i++) {
                x0 = PAD_L + (cum/total)*PW
                cum += len[i]
                x1 = PAD_L + (cum/total)*PW
                y  = PAD_T + PH - (len[i]/maxlen)*PH
                path = path sprintf(" L %.2f,%.2f L %.2f,%.2f", x0, y, x1, y)
                if (n50_x < 0 && cum >= total*0.5) { n50_x=(x0+x1)/2; n50_y=y; n50_len=len[i] }
            }
            path = path sprintf(" L %.2f,%.2f Z", PAD_L+PW, PAD_T+PH)
            printf "<path d=\"%s\" fill=\"#fdf1e0\" stroke=\"#e99e44\" stroke-width=\"1.5\"/>\n", path
            for (p=0; p<=100; p+=25) {
                tx = PAD_L + (p/100)*PW
                printf "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#ddd\" stroke-dasharray=\"2,2\"/>\n", tx, PAD_T, tx, PAD_T+PH
                printf "<text x=\"%.2f\" y=\"%d\" font-size=\"9\" fill=\"#666\" text-anchor=\"middle\">%d%%</text>\n", tx, H-6, p
            }
            for (q=0; q<=4; q++) {
                qy = PAD_T + PH - (q/4)*PH
                qv = int((q/4)*maxlen)
                if (q>0) printf "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#eee\"/>\n", PAD_L, qy, PAD_L+PW, qy
                printf "<text x=\"%d\" y=\"%.2f\" font-size=\"9\" fill=\"#666\" text-anchor=\"end\">%s</text>\n", PAD_L-4, qy+3, qv
            }
            printf "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#999\"/>\n", PAD_L, PAD_T+PH, PAD_L+PW, PAD_T+PH
            printf "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#999\"/>\n", PAD_L, PAD_T, PAD_L, PAD_T+PH
            if (n50_x > 0) {
                printf "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#c0392b\" stroke-width=\"1.2\" stroke-dasharray=\"3,2\"/>\n", n50_x, PAD_T, n50_x, PAD_T+PH
                printf "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3\" fill=\"#c0392b\"><title>N50: %s bp</title></circle>\n", n50_x, n50_y, n50_len
                label_anchor = (n50_x > W - PAD_R - 30) ? "end" : "start"
                label_x = (label_anchor == "end") ? n50_x - 6 : n50_x + 6
                printf "<rect x=\"%.2f\" y=\"%.2f\" width=\"30\" height=\"12\" fill=\"#fff\" opacity=\"0.8\"/>\n", (label_anchor=="end" ? label_x-30 : label_x), PAD_T+2
                printf "<text x=\"%.2f\" y=\"%.2f\" font-size=\"9\" fill=\"#c0392b\" text-anchor=\"%s\">N50</text>\n", label_x, PAD_T+11, label_anchor
            }
        }'
}

_svg_contig_bars() {
    local lengths="$1" cap=30
    [[ -n "${lengths}" ]] || return 0
    local max_len; max_len=$(sort -rn <<< "${lengths}" | head -1)
    [[ "${max_len}" -gt 0 ]] || return 0
    local shown; shown=$(wc -l <<< "${lengths}")
    (( shown > cap )) && shown=${cap}
    sort -rn <<< "${lengths}" | head -"${shown}" | LC_NUMERIC=C awk -v max="${max_len}" -v n="${shown}" '
        BEGIN { W = 560; H = 110; gap = 2; bw = (W-(n-1)*gap)/n }
        {
            i++
            h = ($1/max)*H; if (h < 1) h = 1
            x = (i-1)*(bw+gap); y = H-h
            printf "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#e99e44\"><title>%s bp</title></rect>\n", x, y, bw, h, $1
        }'
}

generate_report() {
    log "Generating final report..."

    local n_contigs total_len n50 largest
    n_contigs=$(grep -c '^>' "${ASSEMBLY}" 2>/dev/null || echo "N/A")

    if command -v assembly-stats &>/dev/null; then
        local _stats
        _stats=$(assembly-stats "${ASSEMBLY}" 2>/dev/null)
        total_len=$(awk '/^sum[ =]|sum[_ ]len|total[_ ]len/{for(i=1;i<=NF;i++){v=$i; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}' <<< "${_stats}")
        n50=$(      awk '/^N50[ =]/{for(i=1;i<=NF;i++){v=$i; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}'         <<< "${_stats}")
        largest=$(  awk '{for(i=1;i<=NF;i++){if($i~/largest/){for(j=i;j<=NF;j++){v=$j; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}}}' <<< "${_stats}")
        total_len="${total_len:-N/A}"
        n50="${n50:-N/A}"
        largest="${largest:-N/A}"
    else
        total_len="N/A (install assembly-stats)"
        n50="N/A"; largest="N/A"
    fi

    local mlst_st="N/A" amr_gene_count="N/A" abricate_gene_count="N/A" pf_hit_count="N/A"
    if ${MLST_DETECTION} && [[ -s "${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv" ]]; then
        mlst_st=$(awk -F'\t' 'NR==1{print $3}' "${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv")
    fi
    if ${AMR_DETECTION} && [[ -s "${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv" ]]; then
        amr_gene_count=$(( $(wc -l < "${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv") - 1 ))
    fi
    if ${ABRICATE_DETECTION} && [[ -s "${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv" ]]; then
        abricate_gene_count=$(( $(wc -l < "${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv") - 1 ))
    fi
    if ${PLASMID_DETECTION}; then
        local pf_dir="${OUTDIR}/reports/plasmidfinder"
        local pf_tab
        pf_tab=$(find "${pf_dir}" -maxdepth 1 -iname 'results_tab.tsv' 2>/dev/null | head -1)
        if [[ -n "${pf_tab}" && -s "${pf_tab}" ]]; then
            pf_hit_count=$(( $(wc -l < "${pf_tab}") - 1 ))
        elif [[ -s "${pf_dir}/data.json" ]]; then
            if command -v python3 &>/dev/null; then
                pf_hit_count=$(python3 -c "
import json
d = json.load(open('${pf_dir}/data.json'))['plasmidfinder']['results']
print(sum(len(s) for db in d.values() for s in db.values() if isinstance(s, dict)))
" 2>/dev/null) || pf_hit_count="N/A (parse error)"
            else
                pf_hit_count="N/A (python3 not found)"
            fi
        fi
    fi

    local quast_tsv="${OUTDIR}/reports/quast/report.tsv"
    local quast_gc="N/A" quast_ns="N/A" quast_rrna="N/A" quast_genes="N/A"
    if [[ -s "${quast_tsv}" ]]; then
        quast_gc=$(awk    -F'\t' '$1=="GC (%)"{print $2}'                                             "${quast_tsv}")
        quast_ns=$(awk     -F'\t' 'index($1,"# N")==1 && index($1,"kbp")>0{print $2; exit}'             "${quast_tsv}")
        quast_rrna=$(awk   -F'\t' 'index($1,"# predicted rRNA genes")==1{print $2; exit}'               "${quast_tsv}")
        quast_genes=$(awk  -F'\t' 'index($1,"# predicted genes")==1{print $2; exit}'                    "${quast_tsv}")
        quast_gc="${quast_gc:-N/A}"; quast_ns="${quast_ns:-N/A}"
        quast_rrna="${quast_rrna:-N/A}"; quast_genes="${quast_genes:-N/A}"
    fi

    local busco_summary busco_lineage="N/A" busco_complete="N/A" busco_single="N/A" busco_dup="N/A" busco_frag="N/A" busco_missing="N/A"
    busco_summary=$(find "${OUTDIR}/reports/busco" -maxdepth 2 -iname 'short_summary*.txt' 2>/dev/null | head -1)
    if [[ -n "${busco_summary}" ]]; then
        busco_lineage=$(sed -n 's/.*lineage dataset is: \([^ ]*\).*/\1/p' "${busco_summary}" | head -1)
        read -r busco_complete busco_single busco_dup busco_frag busco_missing < <(
            sed -n 's/.*C:\([0-9.]*\)%\[S:\([0-9.]*\)%,D:\([0-9.]*\)%\],F:\([0-9.]*\)%,M:\([0-9.]*\)%,n:[0-9]*.*/\1\t\2\t\3\t\4\t\5/p' "${busco_summary}" | head -1
        ) || true
        busco_lineage="${busco_lineage:-N/A}"; busco_complete="${busco_complete:-N/A}"
        busco_single="${busco_single:-N/A}"; busco_dup="${busco_dup:-N/A}"
        busco_frag="${busco_frag:-N/A}"; busco_missing="${busco_missing:-N/A}"
    fi

    local busco_chart_svg
    busco_chart_svg=$(_svg_busco_bar "${busco_single}" "${busco_dup}" "${busco_frag}" "${busco_missing}")
    local busco_chart_html=""
    if [[ -n "${busco_chart_svg}" ]]; then
        busco_chart_html="<h3>BUSCO completeness</h3>
<svg viewBox=\"0 0 560 36\" width=\"100%\" height=\"36\" xmlns=\"http://www.w3.org/2000/svg\">${busco_chart_svg}</svg>
<p style=\"font-size:.85em;color:#666;\"><span style=\"color:#e99e44;\">■</span> single-copy &nbsp;<span style=\"color:#3492c7;\">■</span> duplicated &nbsp;<span style=\"color:#f0c419;\">■</span> fragmented &nbsp;<span style=\"color:#c0392b;\">■</span> missing</p>"
    fi

    local logo_img_html=""
    if [[ -s "${SCRIPT_DIR}/logo.png" ]]; then
        logo_img_html="<img src=\"data:image/png;base64,$(base64 -w0 "${SCRIPT_DIR}/logo.png")\" alt=\"TakiLine logo\" width=\"48\" height=\"48\" style=\"vertical-align:middle;margin-right:.5rem;\">"
    fi

    local contig_lengths
    contig_lengths=$(awk '/^>/{if(seqlen>0) print seqlen; seqlen=0; next}{seqlen+=length($0)}END{if(seqlen>0) print seqlen}' "${ASSEMBLY}" 2>/dev/null || true)
    local contig_chart_svg="" contig_bars_svg=""
    if [[ -n "${contig_lengths}" ]]; then
        contig_chart_svg=$(_svg_contig_chart "${contig_lengths}")
        contig_bars_svg=$(_svg_contig_bars "${contig_lengths}")
    fi
    local contig_chart_html=""
    if [[ -n "${contig_chart_svg}" ]]; then
        contig_chart_html="<h3>Contig length distribution</h3>
<svg viewBox=\"0 0 700 260\" width=\"100%\" height=\"260\" xmlns=\"http://www.w3.org/2000/svg\">${contig_chart_svg}</svg>
<p style=\"font-size:.85em;color:#666;\">Nx plot — contig length (y) at each point of the assembly, ordered largest to smallest (x, % of total length). The steeper the drop, the fewer contigs carry most of the genome; the marked point is N50.</p>"
        if [[ -n "${contig_bars_svg}" ]]; then
            local _n_contigs_for_chart; _n_contigs_for_chart=$(wc -l <<< "${contig_lengths}")
            local _bars_note=""
            (( _n_contigs_for_chart > 30 )) && _bars_note=" (top 30 of ${_n_contigs_for_chart} shown)"
            contig_chart_html+="
<h3>Individual contigs${_bars_note}</h3>
<svg viewBox=\"0 0 560 110\" width=\"100%\" height=\"110\" xmlns=\"http://www.w3.org/2000/svg\">${contig_bars_svg}</svg>
<p style=\"font-size:.85em;color:#666;\">Sorted by length, largest first.</p>"
        fi
    fi

    local checkm_tsv="" checkm_tool="N/A" checkm_completeness="N/A" checkm_contamination="N/A"
    if [[ -s "${OUTDIR}/reports/checkm2/quality_report.tsv" ]]; then
        checkm_tsv="${OUTDIR}/reports/checkm2/quality_report.tsv"; checkm_tool="CheckM2"
    elif [[ -s "${OUTDIR}/reports/checkm/checkm_results.tsv" ]]; then
        checkm_tsv="${OUTDIR}/reports/checkm/checkm_results.tsv"; checkm_tool="CheckM (legacy)"
    fi
    if [[ -n "${checkm_tsv}" ]]; then
        read -r checkm_completeness checkm_contamination < <(
            awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Completeness")c=i; if($i=="Contamination")k=i}} NR==2{print $c"\t"$k}' "${checkm_tsv}"
        ) || true
        checkm_completeness="${checkm_completeness:-N/A}"; checkm_contamination="${checkm_contamination:-N/A}"
    fi

    local kraken_top_pct="N/A"
    [[ -f "${OUTDIR}/logs/.species_pct" ]] && kraken_top_pct=$(cat "${OUTDIR}/logs/.species_pct")

    local busco_flag checkm_completeness_flag checkm_contamination_flag kraken_flag
    busco_flag=$(_flag "${busco_complete}" "<95" true)
    checkm_completeness_flag=$(_flag "${checkm_completeness}" "<90" true)
    checkm_contamination_flag=$(_flag "${checkm_contamination}" ">5" true)
    # kraken_top_pct is legitimately N/A when -S is given without --verify-species (Kraken2 skipped by design).
    kraken_flag=$(_flag "${kraken_top_pct}" "<70")

    # Unpolished ONT indels depress BUSCO/CheckM2 scores independent of true quality.
    local polish_caveat_md="" polish_caveat_html=""
    if [[ "${READ_TYPE}" == "long" && "${HIFI_MODE}" == false && -f "${OUTDIR}/logs/.medaka_unpolished" ]]; then
        polish_caveat_md="> ⚠️ **Unpolished assembly:** Medaka did not produce a consensus (see \`${OUTDIR}/logs/medaka.log\`). Low BUSCO/CheckM2 completeness or elevated contamination below may reflect unpolished ONT indel errors rather than a genuinely bad assembly.
"
        polish_caveat_html="<p style=\"color:#a15c00;\">⚠️ <strong>Unpolished assembly:</strong> Medaka did not produce a consensus (see <code>${OUTDIR}/logs/medaka.log</code>). Low BUSCO/CheckM2 completeness or elevated contamination below may reflect unpolished ONT indel errors rather than a genuinely bad assembly.</p>"
    fi

    local il_reads_before="N/A" il_reads_after="N/A" il_bases_after="N/A" il_q20_after="N/A" il_q30_after="N/A" il_pct_passed="N/A"
    local il_json="${OUTDIR}/qc/${SAMPLE}_fastp.json"
    if [[ -s "${il_json}" ]]; then
        read -r il_reads_before il_reads_after il_bases_after il_q20_after il_q30_after < <(
            python3 -c "
import json
d = json.load(open('${il_json}'))
b = d['summary']['before_filtering']; a = d['summary']['after_filtering']
print(b['total_reads'], a['total_reads'], a['total_bases'], round(a['q20_rate']*100,1), round(a['q30_rate']*100,1), sep='\t')
" 2>/dev/null
        ) || true
        il_reads_before="${il_reads_before:-N/A}"; il_reads_after="${il_reads_after:-N/A}"
        il_bases_after="${il_bases_after:-N/A}"; il_q20_after="${il_q20_after:-N/A}"; il_q30_after="${il_q30_after:-N/A}"
        if [[ "${il_reads_before}" != "N/A" && "${il_reads_after}" != "N/A" ]]; then
            il_pct_passed=$(LC_NUMERIC=C awk -v a="${il_reads_after}" -v b="${il_reads_before}" 'BEGIN{if(b>0) printf "%.1f", (a/b)*100; else print "N/A"}')
        fi
    fi

    local lr_num_reads="N/A" lr_total_bases="N/A" lr_mean_len="N/A" lr_mean_qual="N/A" lr_n50="N/A"
    local nanostats="${OUTDIR}/qc/nanoplot/NanoStats.txt"
    if [[ -s "${nanostats}" ]]; then
        lr_num_reads=$(_nanostat "Number of reads" "${nanostats}");    lr_num_reads="${lr_num_reads:-N/A}"
        lr_total_bases=$(_nanostat "Total bases" "${nanostats}");      lr_total_bases="${lr_total_bases:-N/A}"
        lr_mean_len=$(_nanostat "Mean read length" "${nanostats}");    lr_mean_len="${lr_mean_len:-N/A}"
        lr_mean_qual=$(_nanostat "Mean read quality" "${nanostats}");  lr_mean_qual="${lr_mean_qual:-N/A}"
        lr_n50=$(_nanostat "Read length N50" "${nanostats}");          lr_n50="${lr_n50:-N/A}"
    fi

    # READ_TYPE=long is the only mode where Kraken2 classifies just the long-read pool, so it's the only one where low depth makes the % noisy.
    local kraken_confidence_note_md="" kraken_confidence_note_html=""
    if [[ "${READ_TYPE}" == "long" && "${lr_num_reads}" =~ ^[0-9]+$ ]] && (( lr_num_reads < 5000 )); then
        kraken_confidence_note_md=" (low confidence: only ${lr_num_reads} long reads classified)"
        kraken_confidence_note_html=" <span style=\"color:#a15c00;\">(low confidence: only ${lr_num_reads} long reads classified)</span>"
    fi

    local est_coverage="N/A" coverage_flag="" coverage_display="N/A"
    # Mirrors run_qc()'s _target_cov: trycycler pools ~100x for independent subsamples,
    # every other assembler (incl. hybrid/Illumina) targets the general 20x.
    local coverage_threshold=20
    [[ "${READ_TYPE}" == "long" && "${ASSEMBLER}" == "trycycler" ]] && coverage_threshold=100
    local _gs_bases; _gs_bases=$(genome_size_to_bases "${GENOME_SIZE}" 2>/dev/null) || true
    if [[ -n "${_gs_bases}" ]]; then
        local _bases_used=0
        [[ "${il_bases_after}" != "N/A" ]] && _bases_used=$(( _bases_used + ${il_bases_after%.*} ))
        # Use what filtlong actually kept (fed to the assembler), not NanoPlot's raw pre-filter count.
        local _lr_bases_for_cov="${lr_total_bases}"
        [[ -s "${OUTDIR}/logs/.filtered_long_bases" ]] && _lr_bases_for_cov=$(cat "${OUTDIR}/logs/.filtered_long_bases")
        [[ "${_lr_bases_for_cov}" != "N/A" ]] && _bases_used=$(( _bases_used + ${_lr_bases_for_cov%.*} ))
        if (( _bases_used > 0 )); then
            est_coverage=$(LC_NUMERIC=C awk -v b="${_bases_used}" -v g="${_gs_bases}" 'BEGIN{printf "%.1f", b/g}')
            coverage_flag=$(_flag "${est_coverage}" "<${coverage_threshold}")
            coverage_display="${est_coverage}x (min ${coverage_threshold}x)"
        fi
    fi

    local overall_verdict_emoji="✅" overall_verdict_text="PASS" _f
    for _f in "${busco_flag}" "${checkm_completeness_flag}" "${checkm_contamination_flag}" "${kraken_flag}" "${coverage_flag}"; do
        if [[ "${_f}" == "⚠️ " ]]; then
            overall_verdict_emoji="⚠️"; overall_verdict_text="Review recommended (see flagged metrics below)"
            break
        fi
    done

    local read_qc_section="" read_qc_section_html=""
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then
        read_qc_section+="### Illumina (fastp)

| Metric | Value |
|:---|:---|
| Reads before / after filtering | ${il_reads_before} / ${il_reads_after} |
| Reads passing filter | $(_pct "${il_pct_passed}") |
| Bases after filtering | ${il_bases_after} |
| Q20 / Q30 rate (after filtering) | $(_pct "${il_q20_after}") / $(_pct "${il_q30_after}") |

"
        read_qc_section_html+="<h3>Illumina (fastp)</h3>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Reads before / after filtering</td><td>${il_reads_before} / ${il_reads_after}</td></tr>
<tr><td>Reads passing filter</td><td>$(_pct "${il_pct_passed}")</td></tr>
<tr><td>Bases after filtering</td><td>${il_bases_after}</td></tr>
<tr><td>Q20 / Q30 rate (after filtering)</td><td>$(_pct "${il_q20_after}") / $(_pct "${il_q30_after}")</td></tr>
</table>
"
    fi
    if [[ "${READ_TYPE}" == "long" || "${READ_TYPE}" == "hybrid" ]]; then
        read_qc_section+="### Long reads (NanoPlot)

| Metric | Value |
|:---|:---|
| Reads | ${lr_num_reads} |
| Total bases | ${lr_total_bases} |
| Mean read length | ${lr_mean_len} |
| Mean read quality | ${lr_mean_qual} |
| Read length N50 | ${lr_n50} |

"
        read_qc_section_html+="<h3>Long reads (NanoPlot)</h3>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Reads</td><td>${lr_num_reads}</td></tr>
<tr><td>Total bases</td><td>${lr_total_bases}</td></tr>
<tr><td>Mean read length</td><td>${lr_mean_len}</td></tr>
<tr><td>Mean read quality</td><td>${lr_mean_qual}</td></tr>
<tr><td>Read length N50</td><td>${lr_n50}</td></tr>
</table>
"
    fi
    read_qc_section+="**Estimated input coverage:** ${coverage_flag}${coverage_display}
"
    read_qc_section_html+="<p><strong>Estimated input coverage:</strong> ${coverage_flag}${coverage_display}</p>
"

    local species_report_note="" species_source_val=""
    [[ -f "${OUTDIR}/logs/.species_source" ]] && species_source_val=$(cat "${OUTDIR}/logs/.species_source")
    case "${species_source_val}" in
        kraken2)
            species_report_note=" (Kraken2-detected — review \`${OUTDIR}/reports/kraken2/${SAMPLE}_kraken2_report.tsv\` if unexpected)" ;;
        "user (Kraken2-verified)")
            species_report_note=" (user-informed, Kraken2-verified)" ;;
    esac

    local typing_section="" typing_section_html=""
    if ${MLST_DETECTION} || ${AMR_DETECTION} || ${ABRICATE_DETECTION} || ${PLASMID_DETECTION}; then
        typing_section+="## Typing & Resistance"$'\n\n'
        typing_section_html="<h2>Typing &amp; Resistance</h2>"$'\n'"<ul>"$'\n'
        if ${MLST_DETECTION}; then
            typing_section+="- MLST sequence type: **${mlst_st}** (\`${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv\`)"$'\n'
            typing_section_html+="<li>MLST sequence type: <strong>${mlst_st}</strong> (<code>${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv</code>)</li>"$'\n'
        fi
        if ${AMR_DETECTION}; then
            typing_section+="- AMRFinderPlus: **${amr_gene_count}** gene(s)/mutation(s) (\`${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv\`)"$'\n'
            typing_section_html+="<li>AMRFinderPlus: <strong>${amr_gene_count}</strong> gene(s)/mutation(s) (<code>${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv</code>)</li>"$'\n'
        fi
        if ${ABRICATE_DETECTION}; then
            typing_section+="- abricate (${ABRICATE_DB}): **${abricate_gene_count}** hit(s) (\`${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv\`)"$'\n'
            typing_section_html+="<li>abricate (${ABRICATE_DB}): <strong>${abricate_gene_count}</strong> hit(s) (<code>${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv</code>)</li>"$'\n'
        fi
        if ${PLASMID_DETECTION}; then
            typing_section+="- PlasmidFinder: **${pf_hit_count}** replicon hit(s) (\`${OUTDIR}/reports/plasmidfinder/\`)"$'\n'
            typing_section_html+="<li>PlasmidFinder: <strong>${pf_hit_count}</strong> replicon hit(s) (<code>${OUTDIR}/reports/plasmidfinder/</code>)</li>"$'\n'
        fi
        typing_section+=$'\n'
        typing_section_html+="</ul>"$'\n'
    fi

    local input_files_html=""
    [[ -n "${R1:-}"   ]] && input_files_html+="<li>Illumina R1: <code>${R1}</code></li>"$'\n'
    [[ -n "${R2:-}"   ]] && input_files_html+="<li>Illumina R2: <code>${R2}</code></li>"$'\n'
    [[ -n "${LONG:-}" ]] && input_files_html+="<li>Long reads: <code>${LONG}</code></li>"$'\n'

    cat > "${OUTDIR}/reports/SUMMARY.md" <<EOF

**Overall: ${overall_verdict_emoji} ${overall_verdict_text}**

| Field           | Value                       |
|:----------------|:----------------------------|
| Sample          | ${SAMPLE}                   |
| Species         | ${SPECIES:-not determined}${species_report_note} |
| Read type       | ${READ_TYPE}                |
| Assembler       | ${ASSEMBLER}                |
| HiFi mode       | ${HIFI_MODE}                |

$(  [[ -n "${R1:-}"   ]] && echo "- Illumina R1: \`${R1}\`")
$(  [[ -n "${R2:-}"   ]] && echo "- Illumina R2: \`${R2}\`")
$(  [[ -n "${LONG:-}" ]] && echo "- Long reads:  \`${LONG}\`")


${read_qc_section}

| Metric          | Value        |
|:----------------|:-------------|
| Contigs         | ${n_contigs} |
| Total length    | ${total_len} |
| Largest contig  | ${largest}   |
| N50             | ${n50}       |
| Min contig len  | ${MIN_CONTIG_LENGTH} bp |


| Metric                        | Value |
|:-------------------------------|:------|
| GC content                     | $(_pct "${quast_gc}") |
| N's per 100 kbp                | ${quast_ns} |
| Predicted genes / rRNA genes    | ${quast_genes} / ${quast_rrna} |
| BUSCO lineage                   | ${busco_lineage} |
| BUSCO complete (S:single, D:dup)| ${busco_flag}$(_pct "${busco_complete}") (S:$(_pct "${busco_single}"), D:$(_pct "${busco_dup}")) |
| BUSCO fragmented / missing      | $(_pct "${busco_frag}") / $(_pct "${busco_missing}") |
| ${checkm_tool} completeness     | ${checkm_completeness_flag}$(_pct "${checkm_completeness}") |
| ${checkm_tool} contamination    | ${checkm_contamination_flag}$(_pct "${checkm_contamination}") |
| Kraken2 top-species reads       | ${kraken_flag}$(_pct "${kraken_top_pct}")${kraken_confidence_note_md} |

${polish_caveat_md}
Full sub-reports: \`${OUTDIR}/reports/quast/report.html\`, \`${OUTDIR}/reports/busco/\`, \`${OUTDIR}/reports/checkm2/\` (or \`checkm/\`), \`${OUTDIR}/reports/kraken2/\`.

${typing_section}## Suggested Next Steps

1. **Annotation**     – \`Bakta\` (INSDC-ready) or \`PGAP\` for NCBI submission
2. **AMR screening**  – built in via \`--amr\`/\`--abricate\` (see above if run); or \`ResFinder\` for a third opinion
3. **MLST typing**    – built in via \`--mlst\` (see above if run)
4. **Phylogenetics**  – \`IQ-TREE2\` or \`FastTree\`
5. **Pan-genome**     – \`Panaroo\` (preferred) or \`Roary\`

## Output Files

| Step            | Path                                        |
|:----------------|:--------------------------------------------|
| Final assembly  | \`${ASSEMBLY}\`                             |
| QC reports      | \`${OUTDIR}/qc/\`                        |
| Quality reports | \`${OUTDIR}/reports/\`                   |
| Logs            | \`${OUTDIR}/logs/\`                         |
EOF

    cat > "${OUTDIR}/reports/SUMMARY.html" <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TakiLine Report — ${SAMPLE}</title>
<style>
  body { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; max-width: 900px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; line-height: 1.5; }
  h1 { border-bottom: 3px solid #e99e44; padding-bottom: .3rem; }
  h2 { border-bottom: 1px solid #ddd; padding-bottom: .2rem; margin-top: 2rem; }
  table { border-collapse: collapse; width: 100%; margin: 1rem 0; }
  th, td { border: 1px solid #ddd; padding: .5rem .75rem; text-align: left; }
  th { background: #e99e44; color: #fff; }
  tr:nth-child(even) { background: #f6f6f6; }
  code { background: #eee; padding: .1rem .3rem; border-radius: 3px; font-size: .9em; }
  ul { padding-left: 1.3rem; }
</style>
</head>
<body>
<h1>${logo_img_html}TakiLine Assembly Report</h1>
<div style="border-left:4px solid #e99e44;background:#fff7ec;padding:.6rem 1rem;margin:1rem 0;font-weight:600;">${overall_verdict_emoji} Overall: ${overall_verdict_text}</div>
<table>
<tr><th>Field</th><th>Value</th></tr>
<tr><td>Sample</td><td>${SAMPLE}</td></tr>
<tr><td>Species</td><td>${SPECIES:-not determined}${species_report_note}</td></tr>
<tr><td>Read type</td><td>${READ_TYPE}</td></tr>
<tr><td>Assembler</td><td>${ASSEMBLER}</td></tr>
<tr><td>HiFi mode</td><td>${HIFI_MODE}</td></tr>
</table>

<h2>Input Files</h2>
<ul>
${input_files_html}</ul>

<h2>Read QC</h2>
${read_qc_section_html}

<h2>Assembly Statistics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Contigs</td><td>${n_contigs}</td></tr>
<tr><td>Total length</td><td>${total_len}</td></tr>
<tr><td>Largest contig</td><td>${largest}</td></tr>
<tr><td>N50</td><td>${n50}</td></tr>
<tr><td>Min contig len</td><td>${MIN_CONTIG_LENGTH} bp</td></tr>
</table>

<h2>Assembly QC</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>GC content</td><td>$(_pct "${quast_gc}")</td></tr>
<tr><td>N's per 100 kbp</td><td>${quast_ns}</td></tr>
<tr><td>Predicted genes / rRNA genes</td><td>${quast_genes} / ${quast_rrna}</td></tr>
<tr><td>BUSCO lineage</td><td>${busco_lineage}</td></tr>
<tr><td>BUSCO complete (S:single, D:dup)</td><td>${busco_flag}$(_pct "${busco_complete}") (S:$(_pct "${busco_single}"), D:$(_pct "${busco_dup}"))</td></tr>
<tr><td>BUSCO fragmented / missing</td><td>$(_pct "${busco_frag}") / $(_pct "${busco_missing}")</td></tr>
<tr><td>${checkm_tool} completeness</td><td>${checkm_completeness_flag}$(_pct "${checkm_completeness}")</td></tr>
<tr><td>${checkm_tool} contamination</td><td>${checkm_contamination_flag}$(_pct "${checkm_contamination}")</td></tr>
<tr><td>Kraken2 top-species reads</td><td>${kraken_flag}$(_pct "${kraken_top_pct}")${kraken_confidence_note_html}</td></tr>
</table>
${polish_caveat_html}
<p>Full sub-reports: <code>${OUTDIR}/reports/quast/report.html</code>, <code>${OUTDIR}/reports/busco/</code>, <code>${OUTDIR}/reports/checkm2/</code> (or <code>checkm/</code>), <code>${OUTDIR}/reports/kraken2/</code>.</p>

${busco_chart_html}
${contig_chart_html}

${typing_section_html}
<h2>Suggested Next Steps</h2>
<ol>
<li><strong>Annotation</strong> – <code>Bakta</code> (INSDC-ready) or <code>PGAP</code> for NCBI submission</li>
<li><strong>AMR screening</strong> – built in via <code>--amr</code>/<code>--abricate</code> (see above if run); or <code>ResFinder</code> for a third opinion</li>
<li><strong>MLST typing</strong> – built in via <code>--mlst</code> (see above if run)</li>
<li><strong>Phylogenetics</strong> – <code>IQ-TREE2</code> or <code>FastTree</code></li>
<li><strong>Pan-genome</strong> – <code>Panaroo</code> (preferred) or <code>Roary</code></li>
</ol>
<h2>Output Files</h2>
<table>
<tr><th>Step</th><th>Path</th></tr>
<tr><td>Final assembly</td><td><code>${ASSEMBLY}</code></td></tr>
<tr><td>QC reports</td><td><code>${OUTDIR}/qc/</code></td></tr>
<tr><td>Quality reports</td><td><code>${OUTDIR}/reports/</code></td></tr>
<tr><td>Logs</td><td><code>${OUTDIR}/logs/</code></td></tr>
</table>
</body>
</html>
EOF

    echo -e "\n${GREEN}══════════════════════════════════════${NC}"
    echo -e "${GREEN}       PIPELINE COMPLETE ✓            ${NC}"
    echo -e "${GREEN}══════════════════════════════════════${NC}"
    echo -e "  Sample    : ${SAMPLE}"
    echo -e "  Assembly  : ${ASSEMBLY}"
    echo -e "  Contigs   : ${n_contigs}"
    echo -e "  N50       : ${n50}"
    echo -e "  Report    : ${OUTDIR}/reports/SUMMARY.md (${OUTDIR}/reports/SUMMARY.html)"
    echo -e "${GREEN}══════════════════════════════════════${NC}\n"

    trap - ERR
}

