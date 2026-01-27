# SNPfix

[
[

**SNPfix** — bioinformatics pipeline for transferring SNP coordinates between genome assemblies, including closely related species.

## ✨ Features

- ✅ Transfers SNP coordinates between **different assemblies** and **species**
- ✅ Automatic flanking sequence extraction  
- ✅ BLAST alignment + best hit selection by e-value
- ✅ Chromosome/position sorting & cleaning
- ✅ Handles missing positions (NA)

## 🚀 Quick Start

```bash
git clone https://github.com/afonnikova/snpfix.git
cd snpfix

nextflow run snpfix.nf \
  --fasta_source morex_v1.fa \
  --fasta_target morex_v2.fa \
  --chip barley50k.tsv \
  --flank_size 50
```

## 📋 Input Requirements

| Parameter | Description | Format |
|-----------|-------------|--------|
| `--fasta_source` | Source genome | FASTA + `.fai` |
| `--fasta_target` | Target genome | FASTA + `.fai` |
| `--chip` | SNP coordinates | `Marker<tab>Chr<tab>Position` |
| `--flank_size` | Flanking region | int (default: 50) |

**⚠️ Chromosome names must match exactly** between chip file and genomes.

## 📤 Output Files

```
results/
├── snpfix_final.tsv      # Final SNP coordinates (sorted)
├── flanked.fa           # Flanking sequences  
├── blast_results.tsv    # BLAST hits
├── skipped.log          # Failed markers
└── lengths.tsv          # Flank lengths per marker
```

## 🧪 Example

**Input `barley50k.tsv`:**
```
JHI-Hv50k-2016-94418	2H	422934844
JHI-Hv50k-2016-94431	2H	423044915
```

**Output `snpfix_final.tsv`:**
```
JHI-Hv50k-2016-94418	2H	423012345
JHI-Hv50k-2016-94431	2H	423154678
```

## 💻 Installation

```bash
# Conda (recommended)
conda create -n snpfix -c bioconda nextflow samtools blast pandas
conda activate snpfix

# Docker
docker pull afonnikova/snpfix:latest
```

## ⚙️ Resources

```groovy
params {
    flank_size = 50
    cpus = 4
    memory = '8 GB'
}
```

## 🔬 Tested Species/Arrays

| Species | Array | Assemblies |
|---------|-------|------------|
| Barley | Hv50k | Morex v1→v2 |
| Cattle | BovineSNP50 v3 | UMD3.1→ARS-UCD2 |
| Arabidopsis | 250K | TAIR8→TAIR10 |

## 📊 Performance

- **Size:** 20K (scripts only)
- **Runtime:** 15-40 min (50k SNPs, 4 cpus)  
- **RAM:** 4-16 GB
- **Disk:** 10-15 GB (temp files)


## 📝 Citation

```bibtex
@misc{snpfix2026,
  author = {Afonnikova, Svetlana},
  title  = {SNPfix: SNP coordinate transfer pipeline},
  year   = {2026},
  url    = {https://github.com/afonnikova/snpfix}
}
```

***

**Compact (20K), portable, production-ready.** 🚀