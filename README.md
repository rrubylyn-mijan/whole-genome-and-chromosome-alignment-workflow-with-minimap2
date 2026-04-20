# Whole Genome and Chromosome Alignment Workflow with Minimap2

This workflow describes how to:
1. Split genomes by **chunks** (100 Mb) or by **chromosomes**.
2. Align **chunks** or **chromosomes** to a whole reference genome using **Minimap2**.
3. Combine `.paf` outputs into continuous files.
4. Generate dot plots in **R** for whole-genome and chromosome-level visualizations.
5. Identify primary and non-primary alignments.

---

## 1. Split Genome into Chunks (100 Mb)

```bash
ml samtools/1.19.2

FASTA_DIR=/directory/this/saved/split_wheat_by_chr
CHUNK_SIZE=100000000

mkdir -p /directory/this/saved/split_wheat_by_chunks

for CHROM in $FASTA_DIR/*.fa; do
  CHR_NAME=$(basename $CHROM .fa)
  samtools faidx $CHROM
  CHR_LEN=$(awk '{print $2}' ${CHROM}.fai)
  for START in $(seq 1 $CHUNK_SIZE $CHR_LEN); do
    END=$((START + CHUNK_SIZE - 1))
    if [ $END -gt $CHR_LEN ]; then END=$CHR_LEN; fi
    samtools faidx $CHROM $CHR_NAME:$START-$END > /directory/this/saved/split_wheat_by_chunks/${CHR_NAME}_${START}_${END}.fa
  done
done
```

## 2. Split Genome by Chromosome
```bash
ml samtools/1.19.2
samtools faidx /directory/this/saved/wheat.fasta
mkdir -p /directory/this/saved/split_wheat_by_chr

while read chr; do
  samtools faidx /directory/this/saved/wheat.fasta "$chr" > /directory/this/saved/split_wheat_by_chr/"$chr".fasta
done < chromosomes-wheat.txt
```

## 3. Index Reference with Minimap2
```bash
ml minimap2/2.26
minimap2 -d wheat.mmi wheat.fa
```

## 4. Chunk vs Whole Genome Alignment
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p atlas
#SBATCH --mem=367GB
#SBATCH -J mmi_chunk_wheat
#SBATCH -A genolabswheatphg

ml minimap2/2.26

chunk_dir="/directory/this/saved/split_wheat_by_chunks"
output_dir="/directory/this/saved/minimap2_mmi"
mkdir -p "$output_dir"

for chunk in "$chunk_dir"/*.fasta; do
  base_name=$(basename "$chunk" .fasta)
  minimap2 -cx asm5 /directory/this/saved/wheat.mmi "$chunk" > "$output_dir/${base_name}.paf"
done
```

## 5. Chromosome vs Whole Genome Alignment
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p atlas
#SBATCH --mem=367GB
#SBATCH -J mmi_wheat_chr3A
#SBATCH -A genolabswheatphg

ml minimap2/2.26

chunk_dir="/directory/this/saved/split_wheat_by_chunks"
output_dir="/directory/this/saved/minimap2_mmi_wheat_per_chromosome"
mkdir -p "$output_dir"

for chunk in "$chunk_dir"/chr3A*.fasta; do
  base_name=$(basename "$chunk" .fasta)
  minimap2 -cx asm5 /directory/this/saved/wheat.mmi "$chunk" > "$output_dir/${base_name}.paf"
done
```

## 6. Combine PAF Files
```bash
# Combine all chunk PAFs into one continuous file
ls chr*_*_*.paf | while read x; do
  chr=$(echo "$x" | cut -f 1 -d '_')
  start=$(echo "$x" | cut -f 2 -d '_')
  paste \
    <(cut -f1-4 "$x" | awk -v CHR="$chr" -v START="$start" '{print CHR,$2+START-100000001,$3+START-1,$4+START-1; }' | sed "s/ /\t/g") \
    <(cut -f5- "$x")
done > frj_wheat_chr_1-7.paf

# Normalize reference chromosome names (Column 6): Chr* -> chr*
awk 'BEGIN{OFS="\t"} {$6 = tolower(substr($6,1,1)) substr($6,2); print}' frj_wheat_chr_1-7.paf \
  | sort -k6,6V -k1,1 -k2,2n > frj_wheat_chr_1-7_modified.paf

# Chromosome 3B only
ls chr3B_*_*.paf | while read x; do
  chr=$(echo "$x" | cut -f 1 -d '_')
  start=$(echo "$x" | cut -f 2 -d '_')
  paste \
    <(cut -f1-4 "$x" | awk -v CHR="$chr" -v START="$start" '{print CHR,$2+START-1,$3+START-1,$4+START-1;}' | sed "s/ /\t/g") \
    <(cut -f5- "$x")
done > frj_wheat_chr3B.paf
```

## 7. Dot Plot (Whole Genome)
```bash
library(ggplot2)

align <- read.table("frj_wheat_chr_1-7_modified.paf", header = FALSE, fill = TRUE)
colnames(align)[c(1,3,4,6,8:12)] <- c("Ta_chr","Ta_start","Ta_stop","Ref_chr","Ref_start","Ref_stop","resMatch","length","qual")

align$Ta_pos  <- (align$Ta_start + (align$Ta_stop  - align$Ta_start )/2) / 1e6
align$Ref_pos <- (align$Ref_start + (align$Ref_stop - align$Ref_start)/2) / 1e6
align$perID   <- align$resMatch / align$length

# Use first 21 Ta chromosomes as level order for Ref
Ta_chrs <- unique(align$Ta_chr)[1:21]
align$Ref_chr <- factor(align$Ref_chr, levels = Ta_chrs)

# Filter
align_use <- subset(align, Ta_chr %in% Ta_chrs & !is.na(Ref_chr) & length > 10000 & qual >= 60)

# Compute cumulative positions (example approach)
Ta_chrs_max <- setNames(aggregate(align_use$Ta_pos, list(align_use$Ta_chr), max), c("chr","max"))
Ta_chrs_max$cum <- cumsum(Ta_chrs_max$max)
align_use$Ta_CumPos <- align_use$Ta_pos + (Ta_chrs_max$cum[match(align_use$Ta_chr, Ta_chrs_max$chr)] - Ta_chrs_max$max[match(align_use$Ta_chr, Ta_chrs_max$chr)])

ref_chrs_max <- setNames(aggregate(align_use$Ref_pos, list(align_use$Ref_chr), max), c("chr","max"))
ref_chrs_max$cum <- cumsum(ref_chrs_max$max)
align_use$Ref_CumPos <- align_use$Ref_pos + (ref_chrs_max$cum[match(align_use$Ref_chr, ref_chrs_max$chr)] - ref_chrs_max$max[match(align_use$Ref_chr, ref_chrs_max$chr)])

ggplot(align_use, aes(x = Ta_CumPos, y = Ref_CumPos)) +
  geom_point(cex = 0.15, aes(col = perID, alpha = length)) +
  scale_colour_gradient(high = "#005660", low = "white") +
  labs(title = "Alignment: Wheat query vs Wheat subject",
       x = "Wheat (query, Mb)", y = "Wheat (subject, Mb)") +
  theme_minimal()

ggsave(filename = "wheat_whole_genome_dotplot.png", width = 10, height = 5, dpi = 600)
```

## 8. Dot Plot (Chromosome-Level, Continuous)
```bash
library(ggplot2)
library(dplyr)

align <- read.table("continuous_combined_wheat_chr3B.paf", header = FALSE, fill = TRUE)
colnames(align)[c(1,3,4,6,8:12)] <- c("Ta_chr","Ta_start","Ta_stop","Ref_chr","Ref_start","Ref_stop","resMatch","length","qual")

align$perID <- align$resMatch / align$length
align_use <- subset(align, Ta_chr == "chr3B" & length > 10000 & qual > 10)

align_use <- align_use %>%
  arrange(Ta_start) %>%
  mutate(
    cum_Ta_pos  = cumsum(c(0, diff(Ta_start))),
    cum_Ref_pos = cumsum(c(0, diff(Ref_start)))
  )

ggplot(align_use, aes(x = cum_Ta_pos/1e6, y = cum_Ref_pos/1e6)) +
  geom_point(cex = 0.15, aes(col = perID)) +
  scale_colour_gradient(high = "#005660", low = "#34A798") +
  labs(x = "Query chr3B (Mb)", y = "Subject chr3B (Mb)") +
  theme_minimal()

ggsave(filename = "chr3B_high_reso.png", width = 10, height = 5, dpi = 600)
```

## 9. File Transfer (HPC → Local)
```bash
scp user.name@atlas-login.hpc.msstate.edu:/directory/this/saved/.../continuous\_combined\_chinese\_vs\_sumai3\_chr3B.paf C:\\Users\\user.name\\Documents\\Dotplot\_R\_script\\
```

## 10. Combining PAF Files \& Filtering Alignments
```bash
# Combine per-chromosome PAFs into one
cat wheat-subject_vs_wheat-query_chr7A*.paf > combined_wheat-subject_vs_wheat-query_chr7A.paf
cat wheat-subject_vs_wheat-query_chr7B*.paf > combined_wheat-subject_vs_wheat-query_chr7B.paf
cat wheat-subject_vs_wheat-query_chr7D*.paf > combined_wheat-subject_vs_wheat-query_chr7D.paf
cat combined_wheat-subject_vs_wheat-query_chr*.paf > combined1_wheat-subject_vs_wheat-query.paf

# Identify and count alignment classes (Minimap2 'tp' tag)
grep 'tp:A:S' combined1_wheat-subject_vs_wheat-query.paf > secondary_alignments.paf
grep 'tp:A:I' combined1_wheat-subject_vs_wheat-query.paf > supplementary_alignments.paf
grep -E 'tp:A:S|tp:A:I' combined1_wheat-subject_vs_wheat-query.paf > non_primary_alignments.paf
grep 'tp:A:P' combined1_wheat-subject_vs_wheat-query.paf > primary_alignments.paf

grep -c 'tp:A:S' combined1_wheat-subject_vs_wheat-query.paf   # Secondary count
grep -c 'tp:A:I' combined1_wheat-subject_vs_wheat-query.paf   # Supplementary count
grep -c 'tp:A:P' combined1_wheat-subject_vs_wheat-query.paf   # Primary count

# Uppercase Chr (if needed) in a PAF (start..end part left intact)
sed -E 's/^(chr)([1-7][ABD]):[0-9]+-[0-9]+/\u\1\2/' non_primary_alignments.paf > non_primary_alignments_fixed.paf
```

Maintainer:

Ruby Mijan


