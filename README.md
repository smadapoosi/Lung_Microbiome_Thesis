# Comparison of 16S rRNA Gene Sequencing and Metagenomic Whole Genome Shotgun Sequencing Methods for Machine Learning on the Asthma Microbiome
My undergraduate honors thesis project in microbiology, using machine learning methods to study the lung microbiome of mild-moderate asthmatics for associations with asthma and asthma-related clinical phenotypes.

## Abstract:
Asthma is an obstructive airway disease characterized by episodic and reversible airway hyperresponsiveness that affects over 300 million worldwide (Dougherty et al., 2007; Jia, Zhang, & Lv, 2013; Loftus & Wise, 2016). In addition to genetic factors and environmental exposures, the composition of the lung microbiome is thought to influence the heterogeneity in clinical presentation, treatment, and prognosis of asthma (Huang & Boushey, 2015). The advent of “Next-generation sequencing” and the development of culture-independent tools for the identification of microbes, paired with long-standing machine learning (ML) and robust statistical methods, has transformed the ability of researchers to examine the relationships between the human microbiome and disease (Huang & Boushey, 2015). We sought to examine whether the lung microbiome differs by asthma status or asthma-related clinical parameters, such as the use of inhaled corticosteroids, as well as whether the granularity of metagenomic whole genome shotgun sequencing (WGS) to the species-level is uniquely informative or if similar results can be gleaned using cheaper high-throughput 16S rRNA gene sequencing methods. The microbiota of study participants was hypothesized to differ by asthma status and asthma-related clinical parameters and the species-level granularity of WGS was expected to be uniquely informative, such that a parallel analysis of both WGS and 16S rRNA gene sequencing results was predicted to be the optimal approach to lung microbiome analysis. 16S rRNA gene sequencing data on 22 and WGS data on 16 mild-moderate asthmatics from the University of Michigan Characterization of Adults for Asthma Microbiome Research prospective observational cohort study (NCT02887911) were analyzed using principal component analysis, permutational multivariate ANOVA, goodness-of-fit tests, Fisher scoring, random forests, and support vector machines. Software utilized included packages in R as well as the command-line tools AWK and sed. All statistical and ML methods demonstrated that the composition of the lung microbiome differed by asthma status as well as asthma-related clinical parameters such as BMI and ICS use among asthmatics. 16S rRNA gene sequencing-defined OTU0033, which mapped to Prevotella, was found to significantly differ in relative abundance by asthma status even after p-value correction (p = 0.04). A recursive feature elimination step in the ML workflow was demonstrated to be necessary for optimal classification accuracy but to be nonessential for microbiota feature ranking precision. However, the sample size examined was not adequately large to optimize ML methods or to demonstrate robust differences in microbiota relative abundances so additional analyses using the same workflow on larger data sets are warranted.

(420 words)

## Table of Contents:

## Installation:

## Usage:

## Contributing:

## References:
Dougherty, D. et al. (2007). National Asthma Education and Prevention Program,
   Third Expert Panel on the Diagnosis and Management of Asthma. Bethesda (MD): National Heart, Lung, and Blood Institute (US).    
   https://www.nhlbi.nih.gov/files/docs/guidelines/asthsumm.pdf <br />
Huang, Y. J., Boushey, H. A. (2015, January). The Microbiome in Asthma. Journal of Allergy
    and Clinical Immunology, 135(1), 25-30. doi:10.1016/j.jaci.2014.11.011 <br />
Jia, C. E., Zhang, H. P. Lv, Y. et al. (2013). The Asthma Control Test and Asthma Control
    Questionnaire for assessing asthma control: Systematic review and meta-analysis. The Journal of Allergy and Clinical Immunology,
    131(3), 696-703. doi:doi.org/10.1016/j.jaci.2012.08.023. <br />
Loftus, P. A., Wise, S. K. (2016) Epidemiology of asthma. Current Opinion in Otolaryngology
    & Head and Neck Surgery, 24(3), 245-249. doi:10.1097/MOO.0000000000000262.


