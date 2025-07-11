import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# --- Start of the logical explanation and calculation ---

# Plan: Analyze the relationship between broad-sense heritability (H^2),
# narrow-sense heritability (h^2), and the predictive power of a polygenic score (PGS).

# --- Definitions ---
# Phenotypic Variance (Vp) = Genetic Variance (Vg) + Environmental Variance (Ve)
# Genetic Variance (Vg) = Additive Variance (Va) + Dominance Variance (Vd) + Epistatic Variance (Vi)
# Broad-sense Heritability (H^2) = Vg / Vp
# Narrow-sense Heritability (h^2) = Va / Vp
# A standard PGS from GWAS explains, at best, a proportion of variance equal to h^2.

# --- Problem Setup ---
# We are given a broad-sense heritability H^2 = 0.5.
H_squared = 0.5

# For demonstration, let's assume total phenotypic variance (Vp) is 100 units.
Vp = 100.0

# From the definition of H^2, we can calculate the total genetic variance (Vg).
# Final Equation: Vg = H^2 * Vp
Vg = H_squared * Vp
print(f"Given H^2 = {H_squared} and assuming Vp = {Vp}, total genetic variance (Vg) is calculated as:")
print(f"Vg = {H_squared} * {Vp} = {Vg}")
print("-" * 40)

# The core principle is that Vg is composed of parts: Vg = Va + Vd + Vi.
# Since Vd and Vi cannot be negative, Va must always be less than or equal to Vg.

print("Analyzing scenarios for the composition of genetic variance (Vg = 50.0):")

# --- Scenario 1: All genetic variance is additive ---
# This means there are no dominance (Vd) or epistatic (Vi) effects.
print("\n--- Scenario 1: No non-additive effects ---")
Vd_1 = 0.0
Vi_1 = 0.0
# Final Equation: Va = Vg - Vd - Vi
Va_1 = Vg - Vd_1 - Vi_1
# The max variance explained by the PGS is h^2.
# Final Equation: h^2 = Va / Vp
h_squared_1 = Va_1 / Vp

print(f"Let Vd = {Vd_1} and Vi = {Vi_1}.")
print(f"Additive variance (Va) = {Vg} - {Vd_1} - {Vi_1} = {Va_1}")
print(f"Max PGS variance explained (h^2) = {Va_1} / {Vp} = {h_squared_1:.2f}")
print(f"Result: PGS variance ({h_squared_1:.2f}) is EQUAL to H^2 ({H_squared}).")


# --- Scenario 2: Some genetic variance is non-additive ---
# Let's assume dominance and epistasis account for 20% of the genetic variance.
print("\n--- Scenario 2: Non-additive effects exist ---")
non_additive_variance_2 = 0.20 * Vg
Vd_2 = non_additive_variance_2 / 2
Vi_2 = non_additive_variance_2 / 2
# Final Equation: Va = Vg - Vd - Vi
Va_2 = Vg - Vd_2 - Vi_2
# Final Equation: h^2 = Va / Vp
h_squared_2 = Va_2 / Vp

print(f"Let Vd = {Vd_2} and Vi = {Vi_2}.")
print(f"Additive variance (Va) = {Vg} - {Vd_2} - {Vi_2} = {Va_2}")
print(f"Max PGS variance explained (h^2) = {Va_2} / {Vp} = {h_squared_2:.2f}")
print(f"Result: PGS variance ({h_squared_2:.2f}) is LESS than H^2 ({H_squared}).")


# --- Conclusion ---
print("\n" + "="*40)
print("CONCLUSION FROM DEMONSTRATION")
print("="*40)
print("The variance a standard PGS can explain is h^2 = Va / Vp.")
print("The total genetic variance is H^2 = (Va + Vd + Vi) / Vp.")
print("\nSince Va <= Va + Vd + Vi, it is always true that (Va/Vp) <= (Vg/Vp).")
print(f"This means h^2 <= H^2. Given H^2 = {H_squared}, h^2 must be <= {H_squared}.")
print("\nTherefore, the polygenic score can NOT explain MORE than 50% of the variance.")
print("Statement A is necessarily true.")
print("Statements B and C are not necessarily true, as either scenario is possible.")
print("Statement D is incorrect as epigenetic effects do not force h^2 < H^2.")

# --- End of the logical explanation and calculation ---

# Restore original stdout
sys.stdout = original_stdout
# Get the content of the string buffer
output = string_buffer.getvalue()

# Print the captured output
print(output)
# Finally, print the answer in the required format
print("<<<A>>>")