import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Main logic of the script ---

# 1. Define SNP heritability (h²_SNP).
# This is the theoretical maximum proportion of phenotypic variance
# that can be explained by all SNPs, assuming a perfect model.
# Let's say for a trait, it's 40%.
snp_heritability = 0.40

# 2. Define real-world limiting factors for a Polygenic Score (PGS).
# The performance of a real PGS is always less than the theoretical maximum
# due to practical limitations.

# Factor A: Finite GWAS sample size leads to imperfect SNP effect estimates.
# Let's model this as capturing only 75% of the true signal.
effect_size_accuracy = 0.75

# Factor B: The PGS model itself has methodological imperfections.
# (e.g., choice of SNPs, handling of correlations (LD)).
# Let's model this as a 90% efficiency.
methodological_efficiency = 0.90

# 3. Calculate the actual predictive ability of the PGS (R-squared).
# It's the theoretical maximum reduced by the real-world loss factors.
pgs_r_squared = snp_heritability * effect_size_accuracy * methodological_efficiency

# 4. Print the comparison.
print(f"SNP Heritability (Theoretical Maximum Variance Explained): {snp_heritability}")
print(f"PGS R-squared (Actual Variance Explained by the model): {pgs_r_squared}")
print("\n# --- Comparison --- #")
is_lower = pgs_r_squared < snp_heritability
print(f"Is the PGS predictive ability ({pgs_r_squared:.3f}) lower than SNP heritability ({snp_heritability:.3f})? {is_lower}")
print("\nThis demonstrates that the variance explained by a polygenic score is lower than the total SNP heritability for the phenotype.")

# --- End of main logic ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user
print(output)