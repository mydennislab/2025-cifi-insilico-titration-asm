import yaml
import os


def get_account_for_jobs(wildcards):
    """Simple alternation based on haplotype"""
    if int(wildcards.hap) == 1:
        return "genome-center-grp"
    else:
        return "mydennisgrp"


YAML = config.get("yaml", "/quobyte/mydennisgrp/projects/vole/data_paths.yaml")
FRAC_LABELS = [str(x) for x in range(10, 101, 10)] # 10,20,...,100
# 5, step of 20
FRAC_LABELS = [str(x) for x in range(20, 101, 20)] # 20,40,60,80,100
# FRAC_LABELS = ["100"]  # for quick testing
# FRAC_LABELS = ["20", "40"]  # for quick testing

def label_to_prop(label: str) -> float: return int(label) / 100.0

with open(YAML) as fh:
    TREE = yaml.safe_load(fh)

# Collect all HindIII libraries (skip the 'hic' node)
HINDIII = []
for species, by_enzyme in TREE.items():
    if "HindIII" not in by_enzyme: continue
    for sample, payload in by_enzyme["HindIII"].items():
        if sample == "hic":  # skip in-house Hi-C R1/R2 node
            continue
        HINDIII.append(dict(
            species=species,
            sample=sample,
            cifi_bam=payload["cifi"]["bam"],
            hifi_bam=payload["hifi"]["bam"],
        ))

def samples_list():
    return [s["sample"] for s in HINDIII]


# Absolute paths to your scripts
# CIFI2PE_SCRIPT = "/quobyte/mydennisgrp/mabuelanin/vole_assembly_project/scripts/cifi2pe_full_length_args.py"
CIFI2PE_SCRIPT = "/quobyte/mydennisgrp/mabuelanin/vole_assembly_project/scripts/optimized_cifi2pe_full_length_args.py"
CALN50_JS      = "/quobyte/mydennisgrp/mabuelanin/vole_assembly_project/scripts/calN50.js"

# ---------------- Targets ----------------
rule all:
    input:
        # Pre-scaffolding stats
        expand("stats/{sample}/{label}/summary.tsv",
               sample=samples_list(), label=FRAC_LABELS),
        # Post-scaffolding (yahs) stats
        expand("stats/{sample}/{label}/yahs_summary.tsv",
               sample=samples_list(), label=FRAC_LABELS)

# ---------------- Core steps ----------------

rule hifi_fasta:
    """HiFi BAM -> FASTA (once per sample)"""
    priority: 100
    input:
        bam=lambda w: next(s for s in HINDIII if s["sample"] == w.sample)["hifi_bam"]
    output:
        fa="hifi/{sample}.hifi.fa"
    threads: 8
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        "samtools fasta -@ {threads} {input.bam} > {output.fa}"


rule cifi_fastq_full:
    """CiFi BAM -> full FASTQ (once per sample), collated (name-grouped), 8 threads"""
    priority: 200
    input:
        bam=lambda w: next(s for s in HINDIII if s["sample"] == w.sample)["cifi_bam"]
    output:
        fq="cifi/{sample}.full.fastq"
    threads: 8
    resources:
        mem_mb=4*1024, runtime=4 * 60, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        # collate -> stdout BAM (-O), then convert to FASTQ
        "samtools collate -O -u -@ {threads} {input.bam} | "
        "samtools fastq   -@ {threads} - > {output.fq}"



rule downsample_cifi:
    """Downsample to 100/75/50/25% (100% = copy full FASTQ)"""
    priority: 300
    input:  "cifi/{sample}.full.fastq"
    output: "cifi/{sample}.{label}.fastq"
    params:
        prop=lambda w: label_to_prop(w.label)
    threads: 16
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        """
        if [ "{wildcards.label}" = "100" ]; then
            cp {input} {output}
        else
            seqtk sample -s100 {input} {params.prop} > {output}
        fi
        """

rule cifi2pe_split:
    """CiFi single FASTQ -> HiC-like PE (R1/R2) with HindIII"""
    priority: 400
    input:  "cifi/{sample}.{label}.fastq"
    output:
        r1="cifi2pe/{sample}.{label}_HiC_R1.fastq",
        r2="cifi2pe/{sample}.{label}_HiC_R2.fastq"
    params:
        out="cifi2pe/{sample}.{label}",
        cutter="HindIII"
    threads: 8
    resources:
        mem_mb=32000, runtime=6 * 60, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        "python {CIFI2PE_SCRIPT} --out {params.out} {input} {params.cutter}"

rule hifiasm_dual_scaf:
    """Assemble with hifiasm --dual-scaf (produces hap1/2 ctg GFAs)"""
    priority: 500
    input:
        r1="cifi2pe/{sample}.{label}_HiC_R1.fastq",
        r2="cifi2pe/{sample}.{label}_HiC_R2.fastq",
        hifi="hifi/{sample}.hifi.fa"
    output:
        hap1_gfa="asm/{sample}/{label}/{sample}.{label}.asm.hic.hap1.p_ctg.gfa",
        hap2_gfa="asm/{sample}/{label}/{sample}.{label}.asm.hic.hap2.p_ctg.gfa"
    params:
        pref="asm/{sample}/{label}/{sample}.{label}.asm"
    threads: 32
    resources:
        mem_mb=300*1024, runtime=24 * 60, slurm_partition="high", slurm_account="genome-center-grp"
    shell:
        "hifiasm --dual-scaf -t {threads} -o {params.pref} "
        "--h1 {input.r1} --h2 {input.r2} {input.hifi}"

# ---------------- QC: GFA -> FASTA -> calN50 -> summary ----------------

rule gfa2fa:
    """Convert GFA contigs to FASTA (hap1 or hap2)"""
    priority: 600
    input:
        gfa="asm/{sample}/{label}/{sample}.{label}.asm.hic.hap{hap}.p_ctg.gfa"
    output:
        fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    threads: 4
    resources:
        mem_mb=32000, runtime=4 * 60, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        "gfatools gfa2fa {input.gfa} > {output.fa}"

rule caln50:
    """Run calN50.js (k8) on each hap FASTA"""
    priority: 700
    input:
        fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        n50="stats/{sample}/{label}/hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"

rule summarize_assembly:
    """
    Parse k8 calN50.js output from hap1/2 and write a tidy TSV with:
    sample, fraction, hap, GS, SZ, NN, N50, L50, AU
    """
    priority: 800
    input:
        hap1="stats/{sample}/{label}/hap1.n50.txt",
        hap2="stats/{sample}/{label}/hap2.n50.txt"
    output:
        tsv="stats/{sample}/{label}/summary.tsv"
    threads: 1
    resources:
        mem_mb=2000, runtime=10, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        r'''
        mkdir -p $(dirname {output.tsv})
        echo -e "sample\tfraction\thap\tGS\tSZ\tNN\tN50\tL50\tAU" > {output.tsv}

        parse_one () {{
          IN="$1"; HAP="$2"; SAMP="{wildcards.sample}"; FRAC="{wildcards.label}";
          GS=$(awk -F'\t' '$1=="GS"{{print $2}}' "$IN")
          SZ=$(awk -F'\t' '$1=="SZ"{{print $2}}' "$IN")
          NN=$(awk -F'\t' '$1=="NN"{{print $2}}' "$IN")
          N50=$(awk -F'\t' '$1=="NL" && $2==50{{print $3}}' "$IN")
          L50=$(awk -F'\t' '$1=="NL" && $2==50{{print $4}}' "$IN")
          AU=$(awk -F'\t' '$1=="AU"{{print $2}}' "$IN")
          echo -e "${{SAMP}}\t${{FRAC}}\t${{HAP}}\t${{GS}}\t${{SZ}}\t${{NN}}\t${{N50}}\t${{L50}}\t${{AU}}"
        }}

        parse_one {input.hap1} hap1 >> {output.tsv}
        parse_one {input.hap2} hap2 >> {output.tsv}
        '''


rule porec_nextflow:
    """Run epi2me-labs/wf-pore-c nextflow pipeline for Hi-C contact mapping"""
    priority: 650
    input:
        cifi_bam=lambda w: next(s for s in HINDIII if s["sample"] == w.sample)["cifi_bam"],
        ref_fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        bed="porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed",
        pairs="porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.pairs.gz",
        mcool="porec/{sample}/{label}/hap{hap}/pairs/{sample}.{label}.hap{hap}.mcool",
        hic="porec/{sample}/{label}/hap{hap}/hi-c/{sample}.{label}.hap{hap}.hic",
        report="porec/{sample}/{label}/hap{hap}/wf-pore-c-report.html"
    params:
        outdir="porec/{sample}/{label}/hap{hap}",
        sample_alias="{sample}.{label}.hap{hap}",
        cutter="HindIII",
        nxf_cache="/quobyte/mydennisgrp/mabuelanin/cache/singularity"
    log:
        nf="porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.nextflow.log",
        time="porec/{sample}/{label}/logs/{sample}.{label}.hap{hap}.time.txt"
    threads: 15
    resources:
        mem_mb=150*1024, runtime=24 * 60, slurm_partition="high", slurm_account=lambda w: "mydennisgrp" if int(w.hap) == 1 else "genome-center-grp"
    shell:
        """
        export NXF_SINGULARITY_CACHEDIR="{params.nxf_cache}"
        export NXF_OPTS='-Xms2g -Xmx4g'
        export NXF_OFFLINE='true'
        export TMPDIR="/scratch/tmp/{wildcards.sample}.{wildcards.label}.hap{wildcards.hap}"
        mkdir -p $TMPDIR
        export NXF_TMP="$TMPDIR"
        export NXF_TEMP="$TMPDIR"
        export NXF_EXECUTOR='local'
        
        # Create and cd to output directory to isolate nextflow run
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        /usr/bin/time -v \
            nextflow run epi2me-labs/wf-pore-c \
                -r v1.3.0 \
                -profile singularity \
                --bam {input.cifi_bam} \
                --ref ../../../../{input.ref_fa} \
                --cutter {params.cutter} \
                --out_dir . \
                --threads {threads} \
                --minimap2_settings '-ax map-hifi' \
                --paired_end \
                --bed --pairs --hi_c --mcool --coverage \
                --sample "{params.sample_alias}" \
                -with-report   pipeline_report.html \
                -with-timeline pipeline_timeline.html \
                -with-trace    pipeline_trace.txt \
                -with-dag      pipeline_dag.svg \
                1> ../../../../{log.nf} 2> ../../../../{log.time}
        """

rule index_fa:
    """Index FASTA files for yahs"""
    priority: 625
    input:
        fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa"
    output:
        fai="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa.fai"
    threads: 1
    resources:
        mem_mb=8 * 1024, runtime=60, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        "samtools faidx {input.fa}"


rule yahs_scaffold:
    """Scaffold haplotype assemblies with yahs"""
    priority: 750
    input:
        asm_fa="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa",
        asm_fai="asm/{sample}/{label}/{sample}.{label}.hap{hap}.fa.fai",
        porec_bed="porec/{sample}/{label}/hap{hap}/bed/{sample}.{label}.hap{hap}.bed"
    output:
        scaffolds="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa",
        agp="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.agp",
        bin="yahs/{sample}/{label}/{sample}.{label}.hap{hap}.bin"
    params:
        prefix="yahs/{sample}/{label}/{sample}.{label}.hap{hap}"
    log:
        "yahs/{sample}/{label}/logs/{sample}.{label}.hap{hap}.yahs.log"
    threads: 8
    resources:
        mem_mb=200*1024, runtime=48 * 60, slurm_partition="high", slurm_account=get_account_for_jobs
    shell:
        """
        yahs \
            -q 0 \
            --no-contig-ec \
            -o {params.prefix} \
            -v 1 \
            {input.asm_fa} \
            {input.porec_bed} \
            2>&1 | tee {log}
        """

rule yahs_caln50:
    """Run calN50.js (k8) on each yahs scaffold FASTA"""
    priority: 850
    input:
        fa="yahs/{sample}/{label}/{sample}.{label}.hap{hap}_scaffolds_final.fa"
    output:
        n50="stats/{sample}/{label}/yahs_hap{hap}.n50.txt"
    threads: 1
    resources:
        mem_mb=4000, runtime=30, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        "mkdir -p $(dirname {output.n50}); "
        "k8 {CALN50_JS} -L2.3g {input.fa} > {output.n50}"


rule summarize_yahs:
    """
    Parse k8 calN50.js output from yahs scaffolds hap1/2 and write a tidy TSV with:
    sample, fraction, hap, GS, SZ, NN, N50, L50, AU
    """
    priority: 900
    input:
        hap1="stats/{sample}/{label}/yahs_hap1.n50.txt",
        hap2="stats/{sample}/{label}/yahs_hap2.n50.txt"
    output:
        tsv="stats/{sample}/{label}/yahs_summary.tsv"
    threads: 1
    resources:
        mem_mb=2000, runtime=10, slurm_partition="high", slurm_account="mydennisgrp"
    shell:
        r'''
        mkdir -p $(dirname {output.tsv})
        echo -e "sample\tfraction\thap\tGS\tSZ\tNN\tN50\tL50\tAU" > {output.tsv}

        parse_one () {{
          IN="$1"; HAP="$2"; SAMP="{wildcards.sample}"; FRAC="{wildcards.label}";
          GS=$(awk -F'\t' '$1=="GS"{{print $2}}' "$IN")
          SZ=$(awk -F'\t' '$1=="SZ"{{print $2}}' "$IN")
          NN=$(awk -F'\t' '$1=="NN"{{print $2}}' "$IN")
          N50=$(awk -F'\t' '$1=="NL" && $2==50{{print $3}}' "$IN")
          L50=$(awk -F'\t' '$1=="NL" && $2==50{{print $4}}' "$IN")
          AU=$(awk -F'\t' '$1=="AU"{{print $2}}' "$IN")
          echo -e "${{SAMP}}\t${{FRAC}}\t${{HAP}}\t${{GS}}\t${{SZ}}\t${{NN}}\t${{N50}}\t${{L50}}\t${{AU}}"
        }}

        parse_one {input.hap1} hap1 >> {output.tsv}
        parse_one {input.hap2} hap2 >> {output.tsv}
        '''