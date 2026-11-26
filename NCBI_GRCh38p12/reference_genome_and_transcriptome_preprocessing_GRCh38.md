# Reference Genome and Transcriptome Preprocessing (GRCh38)

Step-by-step instructions for formatting NCBI genome sequences (FASTA) and gene annotations (GFF3) for downstream analysis.

---

## 1. Reference Genome (FASTA) Preparation

### 1-1. Rename Chromosomes
NCBI genome files use `[chromosome_ID].[version]` identifiers such as `NC_000001.xx`. Convert these to `chr1`, `chr2`, …, `chrX`, `chrY`, and `chrM`.
*(The conversion method is the same as the Perl one-liner used for GFF3 in Section 2-1.)*

### 1-2. Hard-Mask Pseudoautosomal Regions (PAR)
Based on NCBI definitions, hard-mask the PAR regions on `chrY` with `N` using `bedtools maskfasta`.

* **Reference:** [NCBI GRC Human](https://www.ncbi.nlm.nih.gov/grc/human) (PAR coordinates are shared between GRCh38.p12 and GRCh38.p14)

| Name | Chr | Start | Stop |
| :--- | :---: | :---: | :---: |
| PAR#1 | X | 10,001 | 2,781,479 |
| PAR#2 | X | 155,701,383 | 156,030,895 |
| PAR#1 | Y | 10,001 | 2,781,479 |
| PAR#2 | Y | 56,887,903 | 57,217,415 |

```bash
fa=GRCh38.p12_chrxym.fa
bed=par_GRCh38p12.bed

# Mask chrY PAR regions with N (done manually for GRCh38.p12)
bedtools maskfasta -fi ${fa} -bed ${bed} -fo ${fa%.fa}_PARN.fa
```
> **Note:** In the original notes, the output filename used `${fa%.bed}_PARN.fa`; this has been corrected to `${fa%.fa}_PARN.fa` to properly strip the `.fa` extension.

---

## 2. Transcript Reference (GFF3) Preparation

### 2-1. Convert Chromosome Names (`NC_` to `chr`)
Replace `NC_` accession numbers with `chr` names and filter out unplaced/alt contigs.

```bash
gzip -dc ref_GRCh38.p12_top_level.gff3.gz | perl -ne '
  s/NC_000001\.\d+/chr1/;
  s/NC_000002\.\d+/chr2/;
  s/NC_000003\.\d+/chr3/;
  s/NC_000004\.\d+/chr4/;
  s/NC_000005\.\d+/chr5/;
  s/NC_000006\.\d+/chr6/;
  s/NC_000007\.\d+/chr7/;
  s/NC_000008\.\d+/chr8/;
  s/NC_000009\.\d+/chr9/;
  s/NC_000010\.\d+/chr10/;
  s/NC_000011\.\d+/chr11/;
  s/NC_000012\.\d+/chr12/;
  s/NC_000013\.\d+/chr13/;
  s/NC_000014\.\d+/chr14/;
  s/NC_000015\.\d+/chr15/;
  s/NC_000016\.\d+/chr16/;
  s/NC_000017\.\d+/chr17/;
  s/NC_000018\.\d+/chr18/;
  s/NC_000019\.\d+/chr19/;
  s/NC_000020\.\d+/chr20/;
  s/NC_000021\.\d+/chr21/;
  s/NC_000022\.\d+/chr22/;
  s/NC_000023\.\d+/chrX/;
  s/NC_000024\.\d+/chrY/;
  s/NC_012920\.\d+/chrM/;
  print if /^chr/;
' | gzip > GRCh38p12_chrxym.gff3.gz
```

### 2-2. Extract `gene` Records Only
Extract only rows where the 3rd column (feature type) is `gene`.

```bash
gzip -dc GRCh38p12_chrxym.gff3.gz | \
  perl -F'\t' -anle 'print $_ if $F[2] eq "gene";' \
  > GRCh38p12_w_pseudo_v1_chrxym_gene.gff3
```

### 2-3. Append Parent `gene` IDs
Append the parent `gene` ID to the end of all child feature lines (`transcript`, `exon`, etc.) so that all features linked to `pseudo=true` genes can be removed together in Step 2-4.

```bash
perl yab_gff_parse_add_topid_v4b.pl \
  -f=GRCh38p12_chrxym.gff3 \
  -type=topid,type,alldata,pseudo \
  -o=GRCh38p12_chrxym_topid_type.gff3
```

### 2-4. Remove Pseudogenes (`pseudo=true`) and Unwanted Features
Filter out lines containing `pseudo=true` or `partial=true`, as well as `tRNAscan-SE` entries and `rRNA`, `tRNA`, `region`, `match`, and `CDNA_match` features.

```bash
perl -anle '
  print $_ unless $_ =~ /(pseudo=true)|(partial=true)/
    or $F[1] eq "tRNAscan-SE"
    or $F[2] =~ /^(rRNA)|(tRNA)|(region)|(match)|(CDNA_match)$/;
' GRCh38p12_chrxym_topid_type.gff3 > GRCh38p12_wo_pseudo_v1_chrxym_topid_type.gff3
```

### 2-5. Remove Orphaned Records Missing a Parent `gene`
In some `Curated Genomic` entries, the `gene` record was removed while its child `exon` records remained. Trace the `Parent` hierarchy and discard any records that do not resolve to a `gene` ID.

```bash
perl -F'\t' -anle '
  ($id)  = ($F[8] =~ /ID=([^;]+)/)[0];
  ($pid) = ($F[8] =~ /Parent=([^;]+)/)[0];
  if ($id =~ /gene/) {
    $g{$id}++;
    $p{$id} = $id;
  } else {
    $p{$id} = $pid;
    while ($p{$id} ne $p{$p{$id}} and $p{$id}) {
      $p{$id} = $p{$p{$id}};
    }
  }
  if ($p{$id}) {
    print "$_";
  } else {
    print "VOID\t$_";
  }
' GRCh38p12_wo_pseudo_v1_chrxym_topid_type.gff3 | \
  grep -v -e "^VOID\t" > GRCh38p12_wo_pseudo_v1_chrxym_topid_type2.gff3
```

### 2-6. Strip Temporary Tags (`;Top_id` and Beyond)
Remove the temporary `;Top_id` attribute appended in Step 2-3 to restore valid GFF3 formatting.

```bash
perl -anle '$_ =~ s/;Top_id.+$//; print $_;' \
  GRCh38p12_wo_pseudo_v1_chrxym_topid_type2.gff3 \
  > GRCh38p12_wo_pseudo_v1_chrxym.gff3
```

### 2-7. Remove `chrY` Genes Within PAR Regions
Manually delete gene annotations located within the PAR regions on `chrY` to complete the reference preparation.

### 2-8. Convert GFF3 to GTF
Convert the cleaned GFF3 file into GTF format to complete the reference preparation.

```bash
perl yab_gff2gtf_v2.pl \
  -f=GRCh38p12_wo_pseudo_v1_chrxym.gff3 \
  -o=GRCh38p12_wo_pseudo_v1_chrxym.gtf
```