**Transcriptomic insights into developmental and nutritional shifts in Huhu (Prionoplus reticularis) grub larval stages**

*Bennett Henzeler1,6,# , Ruchita Rao Kavle2,7,#, Ngoni Faya4,8, Pascal Giehr5, Alaa El-Din Ahmed Bekhit1, Alan Carne3, Corinna Kersten5, Sabine Schneider5 and Dominic Agyei1,9*

Correspondence:dominic.agyei@monash.edu, bennett.henzeler@cup.lmu.de

The Huhu grub (Prionoplus reticularis; Order: Coleoptera, Family: Cerambycidae, Subfamily: Prioninae - hereafter referred to as P. reticularis) is an edible beetle larva endemic to Aotearoa New Zealand and has long been part of Māori food traditions. Despite its cultural and nutritional relevance, molecular information for this species remains scarce. In this study, we present the first de novo transcriptome assembly for P. reticularis and examine stage‑associated transcriptional differences between small and large feeding larvae using whole‑larva RNA sequencing. Distinct expression patterns were observed between the two larval size classes. Larger larvae showed increased expression of genes involved in protein synthesis, mitochondrial activity, and central metabolic processes, consistent with elevated biosynthetic demand during growth. Enrichment analyses identified over‑representation of functions related to ribosome biogenesis, oxidoreductase activity, and ATP binding, alongside conserved metabolic pathways such as glycolysis, the tricarboxylic acid cycle, and amino‑acid biosynthesis. Differences were also observed in transcripts associated with membrane transport and structural components, suggesting broader physiological shifts accompanying larval growth. Because RNA was extracted from whole larvae, these results represent organism‑level transcriptional differences and do not resolve tissue‑specific regulation. Functional and nutritional interpretations should therefore be viewed as hypothesis‑generating rather than definitive. Nevertheless, this work provides a foundational transcriptomic resource for P. reticularis and offers initial molecular insight into physiological differences within the edible larval phase. The dataset establishes a basis for future tissue‑resolved transcriptomic, proteomic, and metabolomic studies aimed at assessing nutritional properties, bioactive potential, and allergenic considerations of Huhu grub proteins.

**The scripts utilized in this study are provided as supplementary materials to ensure reproducibility.**

## Reusing the shell scripts

Files whose names begin with a number and a space are Bash scripts intended for a
SLURM-managed Linux cluster. They are snapshots of the study workflow, not a portable
workflow manager. Read the header of each script before use: it documents the purpose,
software requirements, expected inputs and outputs, an example invocation, and the
values that normally need adapting.

Submit a script from the directory containing its expected inputs, for example:

```bash
sbatch "2. Trimmomatic"
```

Before submitting on another system:

1. Install or load the named command-line tools and container runtime.
2. Replace sample-specific filenames, relative paths, absolute `/mnt/nfs/...` paths,
   container image locations, and database locations.
3. Adjust `#SBATCH` CPU, memory, time, GPU, and node settings for the local cluster.
4. Confirm that command thread counts do not exceed `--cpus-per-task`.
5. Create a fresh output directory or check the tool's overwrite behavior before a
   rerun, and retain tool/database versions alongside the results.

Some supplementary files contain several analyses in one shell script. SLURM reads
`#SBATCH` directives only before the first executable command, so later directive
blocks do **not** create new jobs or change the active allocation. These cases are
called out in their headers; split them into separate scripts when independent
scheduling or different resources are required.

The numeric filename prefix records the study's broad analysis order. It does not by
itself guarantee that every output path feeds directly into the next numbered file;
the R scripts and manual staging steps are also part of the published analysis.
