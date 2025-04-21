[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15257314.svg)](https://doi.org/10.5281/zenodo.15257314)

# oral-microbiome-periodontal-study

This repository contains the R code and associated input/output files used in the analysis for the paper:

**"The impact of periodontal therapy on the oral microbiome and lung function: an intervention study"**

## Contents

- `script/R_script.R`: Full R pipeline for processing, analysis, and plotting.
- `data/`: Preprocessed OTU table and metadata used in the R pipeline.
- `output/`: Final figures generated for the manuscript and supplementary materials.

## Reproducibility

This repository includes:
- Alpha and beta diversity analysis
- Differential abundance testing (MaAsLin2 and ANCOM-BC)
- Visualizations including boxplots, heatmaps, ordination plots

## Citation

Please cite the paper if using this code or analysis workflow in your research.

---

## Abstract

**Background**: The oral cavity hosts over 700 bacterial species, maintaining a microbial balance crucial for health. Disruptions can lead to periodontitis, which has been linked to systemic conditions, including respiratory diseases. This study examines the effects of periodontal therapy on the subgingival microbiome and airway resistance.

**Methods**: Fifty-seven periodontitis patients (never-smokers, no comorbidities) with high biofilm accumulation and gingival inflammation (plaque score; bleeding on probing >50%) underwent full-mouth periodontal disinfection therapy. Subgingival plaque, pooled from the deepest pockets in each quadrant, and airway resistance (Rrs) were assessed at baseline (T0) and six weeks post-treatment (T1) using the Forced Oscillation Technique. Subgingival plaque samples were analysed with Shotgun sequencing. Alpha- and beta-diversity were evaluated with Wilcoxon rank tests and PERMANOVA, while differentially abundant taxa were identified using ANCOM-BC2. Longitudinal associations between bacterial abundance and airway resistance were assessed using a mixed-effects model (MaAsLin 2.0-package in R). Rrs, age, and BMI was included as fixed effects while accounting for repeated measures.

**Results**: Participants (mean age: 36 years, 58% females) had stage I (21%) or stage II (79%) periodontitis. Individuals with the highest (>median) baseline resistance at 5 Hz (R5), indicating decreased airways patency, had an increased abundance of periodontitis-related taxa, including Prevotella, Porphyromonas, and Tannerella, compared to those with R5 below the median. In contrast, individuals with R5 below the median exhibited a higher abundance of ‘good-periodontal-health’-associated taxa such as Actinomyces and Rothia. At follow-up, Rrs at 11Hz and 19Hz decreased by 4.7% and 5.4%, respectively (p<0.05). Additionally, alpha-diversity significantly declined, and beta-diversity also changed significantly (both p<0.001). The microbiome composition shifted toward health-associated species post-therapy, including Actinomyces oris, A. naeslundii, Rothia dentocariosa, and Lautropia mirabilis. MaAsLin2 identified positive associations between airway resistance and increased abundance of genera such as Prevotella, Tannerella, Treponema, implicated in periodontal inflammation and dysbiosis.

**Conclusion**: Periodontal therapy improves airway function and promotes a healthier microbiome within six weeks. The link between airway resistance and microbial composition highlights the potential role of oral health in respiratory function. These findings emphasize the need for personalized periodontal care, particularly for individuals at risk of respiratory conditions.

