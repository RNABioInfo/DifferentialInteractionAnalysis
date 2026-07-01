# Differential RNA-RNA Interaction Enrichment Analysis

Quarto/R workflow for differential enrichment analysis of RNAnue direct-duplex
interaction counts. Treatment samples are compared with ligation controls;
no-ligation samples are QC-only by default.

## Statistical Question

Primary test:

> Are RNAnue interaction counts higher in treatment samples than in the
> ligation-enabled control background?

This is differential enrichment evidence, not direct proof that an RNA-RNA
interaction is native. Strong calls should also show RNAnue native support,
replicate consistency, clean no-ligation QC, and plausible coverage profiles.

## Sample Roles

Metadata must define sample roles explicitly.

| Role | Use |
| --- | --- |
| `treatment` | Case group for the primary model. |
| `ligation_control` | Ligation-enabled background/control group. |
| `no_ligation_control` | QC-only chemistry/background control. |
| `other_qc` | Additional QC-only samples. |

No-ligation controls are excluded from the primary contrast because including
them would mainly test whether ligation occurred.

## Normalization And Models

The workflow fits negative-binomial count models with edgeR and/or DESeq2.
Counts are normalized with sample exposure and, when available, an
interaction-specific transcript-pair background offset.

```text
log_offset[interaction, sample] =
  log(sample_split_exposure[sample]) +
  log(pair_background_weight[interaction, sample])
```

Key points:

- Sample exposure defaults to RNAnue split-read counts from
  `*_read_counts_summary.tsv`.
- Pair-background offsets use `*_contiguous_transcript_counts.tsv`.
- `pair_background_available` reports whether both interaction arms could use
  the pair-background term in each model sample.
- edgeR keeps numeric RNAnue counts.
- DESeq2 receives rounded counts and writes `DESeq2_rounding_diagnostics.tsv`.
- DESeq2 reports the unshrunken Wald-model LFC as the tested effect size.
  `log2FoldChange_shrunken` is retained separately for display/context.
- Model warnings, including DESeq2/apeglm optimizer warnings, are written to
  `analysis_warnings.tsv`.

Sample QC separates read depth from interaction yield. `mapped_read_count` uses
RNAnue `mapped_reads` when present, otherwise `splits + singletons`.
`interaction_library_sum` is the per-sample sum of the interaction-count matrix.

## Annotation And Target Views

GFF post-processing annotates RNAnue BEDPE arms by coordinate overlap. Partner
class comes from the selected GFF-overlap feature, not from RNAnue arm IDs.

Target-gene views match configured IDs against common GFF attributes (`ID`,
`Name`, `gene`, `locus_tag`, `Alias`, `Parent`) and RNAnue arm feature IDs.
This keeps target lookup usable when RNAnue and the external GFF agree by
coordinates but differ by identifiers.

Important limits:

- Results depend on GFF completeness and coordinate compatibility.
- Broad features can compete with more specific annotations.
- Seqname version mismatches may need normalization.
- Unmatched arms remain `unannotated`; RNAnue IDs are not used as fallback
  partner classes.
- Partner-type matrices summarize candidate classes, not regulatory proof.

## Batch Correction

Known batch covariates can be added as fixed effects:

```yaml
batch_correction:
  mode: auto
  columns: [batch]
  on_confounded: fail
  visualization_remove_batch: true
```

In `auto` mode, a batch column is used only when it is present, complete, has at
least two levels, and gives a full-rank design. Confounded batch columns fail
before model fitting and write `design_diagnostics.tsv` and `batch_balance.tsv`.

Batch correction is model-side only. Optional batch-removed values are used for
PCA visualization, not as edgeR or DESeq2 input.

## Configure Inputs

Copy the example parameter file and edit paths:

```sh
cp params.example.yml params.local.yml
```

Important fields:

- `analysis_method`: `edgeR`, `DESeq2`, or `both`.
- `metadata_file`: CSV with `sample_id`, `condition`, `role`, and `include`.
- `counts_file`: RNAnue `complete_super_interaction_transcript_counts.gct`.
- `interactions_file`: RNAnue `complete_super_interaction_regions.bedpe`.
- `rnanue_results_dir`: RNAnue result root used to discover auxiliary outputs.
- `detect_dir`, `analyze_dir`, `postprocess_dir`: optional explicit overrides.
- `primary_contrast`: roles used in the main model.
- `batch_correction`: known-batch covariate handling.

Bundled Hyderabad dark example roles:

- `01300`, `01187`, `01507`: `treatment`
- `01302`, `01509`: `ligation_control`
- `01190`, `01193`, `01304`: `no_ligation_control`

## Setup

```sh
conda env update -f environment.yml --prune
conda activate differential-interaction-analysis
"${CONDA_PREFIX}/bin/Rscript" code/check_environment.R
```

For the pinned Apple Silicon environment:

```sh
conda create -n differential-interaction-analysis --file environment.osx-arm64.lock.txt
conda activate differential-interaction-analysis
```

## Render

```sh
bash code/quarto.sh render differential_interaction_analysis.qmd --execute-params params.local.yml
```

Preview:

```sh
bash code/quarto.sh preview differential_interaction_analysis.qmd --execute-params params.local.yml
```

## Test

```sh
"${CONDA_PREFIX}/bin/Rscript" tests/run_smoke_tests.R
```

The smoke test uses `tests/fixtures/synthetic_rnanue` and covers role
validation, no-ligation QC, offset construction, edgeR, DESeq2, method
concordance, batch-design validation, post-processing outputs, and deterministic
overwrites.

## Outputs

The report writes deterministic TSV/BED outputs to `out_dir`. Existing files
with the same names are overwritten.

- `edgeR_results.tsv`
- `DESeq2_results.tsv`
- `method_concordance.tsv` when `analysis_method: both`
- `significant_interactions_edgeR.bedpe`
- `significant_interactions_DESeq2.bedpe`
- `high_confidence_interactions.tsv`
- `high_confidence_interactions.bedpe`
- `high_confidence_interaction_regions.bed`
- `high_confidence_interaction_regions.tsv`
- `sample_qc.tsv`
- `analysis_warnings.tsv`
- `normalization_diagnostics.tsv`
- `design_diagnostics.tsv`
- `batch_balance.tsv`
- `candidate_stability.tsv`
- `significant_normalized_counts.tsv`
- `dia_rnanue_padj_comparison.tsv`
- `dia_rnanue_discordance.tsv`
- `rnanue_metric_group_summary.tsv`
- `partner_recurrence.tsv`
- `interaction_arm_gff_annotations.tsv`
- `interaction_partner_type_pairs.tsv`
- `partner_type_confusion.tsv`
- `target_gene_interactions.tsv`
- `target_gene_interactions.bedpe`
- `gff_annotation_diagnostics.tsv`
- `pca_diagnostics.tsv`
- `pca_scores_<set>.tsv`
- `no_ligation_qc.tsv` when no-ligation controls are present
- `DESeq2_rounding_diagnostics.tsv`

## High-Confidence Tier

By default, high-confidence interactions must pass all gates:

- edgeR significant at `padj_thresh` with positive LFC.
- DESeq2 significant at `padj_thresh` with positive LFC.
- RNAnue native minimum `padj_value <= 0.1`.
- no-ligation maximum count equal to `0`.
- pair-background offset available for all model samples.
- coverage profile not exclusively `broad_diffuse`.

Tune these thresholds in the optional `high_confidence` block.

High-confidence outputs:

- `high_confidence_interactions.tsv`: evidence table.
- `high_confidence_interactions.bedpe`: paired interaction rows.
- `high_confidence_interaction_regions.bed`: BED6-style arm regions.
- `high_confidence_interaction_regions.tsv`: arm regions with annotations.

## Interpretation Boundaries

Method-concordant significant interactions are stronger candidates than
method-specific calls, but still need biological review.

Watch for:

- sparse rows driven by one replicate,
- partner RNAs with large abundance shifts,
- high no-ligation signal,
- sample-only offset fallback,
- repeated local clusters sharing an abundant partner.

Useful diagnostics:

- `candidate_stability.tsv`: replicate spread and single-sample dominance.
- `significant_normalized_counts.tsv`: sample-level normalized counts.
- `dia_rnanue_padj_comparison.tsv`: DIA FDR versus RNAnue treatment padj.
- `dia_rnanue_discordance.tsv`: shared, DIA-only, RNAnue-only, and unsupported
  classes.
- `partner_recurrence.tsv`: recurrent partner RNAs.
- PCA and heatmap outputs: sample separation and called-feature structure.

Use the report as a ranked post-processing layer on RNAnue results, not as a
replacement for RNAnue native interaction statistics.

## References

- edgeR: <https://bioconductor.org/packages/release/bioc/html/edgeR.html>
- DESeq2: <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>
- limma `removeBatchEffect`:
  <https://rdrr.io/bioc/limma/man/removeBatchEffect.html>
- diffHic: <https://bioconductor.org/packages/release/bioc/html/diffHic.html>
- TMM normalization:
  <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25>
