# Define hypothetical values for SNP heritability and PGS predictive ability
# SNP heritability (h2_snp) represents the total phenotypic variance
# explained by all considered SNPs. It's a theoretical maximum.
# Let's assume for a given trait, it's 50%.
h2_snp = 0.50

# The predictive ability of a Polygenic Score (PGS R-squared) is the variance
# explained by a predictive model built from GWAS data. Due to finite sample
# sizes in the GWAS, imperfect model assumptions, and other statistical noise,
# the PGS can only capture a fraction of the total SNP heritability.
# A realistic value for a well-powered PGS might be 10-15%.
pgs_r_squared = 0.12

# The core concept is that pgs_r_squared <= h2_snp.
# We can print this relationship.
print(f"The SNP heritability (h2_snp) for a trait is the theoretical maximum proportion of variance explainable by SNPs.")
print(f"A hypothetical h2_snp is: {h2_snp}")
print("\nA Polygenic Score (PGS) is a predictive model based on GWAS data. Its predictive ability (R-squared) is the actual variance it explains in a new sample.")
print(f"A realistic R-squared for the PGS might be: {pgs_r_squared}")
print(f"\nAs such, the PGS predictive ability of {pgs_r_squared} is necessarily lower than the SNP heritability of {h2_snp} because the PGS is an imperfect model trying to capture that total heritable signal.")
