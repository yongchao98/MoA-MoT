import sys
# Redirect print to stderr to show the thinking process, as the final answer will be printed to stdout.
_print = print
def print(s):
    _print(s, file=sys.stderr)

# Plan:
# 1. Define the core concepts and the given value for broad-sense heritability (H^2).
# 2. Establish the key inequality relating the variance explained by a PGS (R^2_PGS),
#    narrow-sense heritability (h^2), and broad-sense heritability (H^2).
# 3. Demonstrate with two scenarios: one with only additive genetic effects, and one with non-additive effects.
# 4. Check which statement holds true for both scenarios.
# 5. Output the final logical conclusion and the final answer.

# Step 1: Define concepts and givens
H_sq = 0.5
print("--- Problem Setup ---")
print(f"Given: Broad-sense heritability (H^2) = {H_sq}")
print("H^2 = (Total Genetic Variance) / (Total Phenotypic Variance)")
print("This is the absolute maximum proportion of variance that can be explained by genetics.")
print("\n--- Key Concepts ---")
print("Narrow-sense heritability (h^2) = (Additive Genetic Variance) / (Total Phenotypic Variance)")
print("A standard Polygenic Score (PGS) primarily captures additive effects.")
print("Therefore, the maximum variance a PGS can explain (R^2_PGS) approaches h^2.")

# Step 2: Establish the key inequality
print("\n--- Core Logical Relationship ---")
print("From the definitions, we have the following inequality:")
print("R^2_PGS <= h^2 <= H^2")
print(f"Substituting the given value for H^2, we get:")
print(f"R^2_PGS <= h^2 <= {H_sq}")
print("This means the variance explained by the PGS cannot be more than 50%.")
print("\n--- Analysis of Answer Choices using Scenarios ---")

# Step 3: Scenario-based analysis
# For simplicity, we assume Total Phenotypic Variance (Vp) = 1.
# So, Vg = H^2 = 0.5 and Va = h^2.

# Scenario 1: All genetic variance is additive.
print("\nSCENARIO 1: No non-additive genetic effects (dominance, epistasis).")
h_sq_1 = 0.5
H_sq_1 = 0.5
max_R_sq_PGS_1 = h_sq_1
print(f"In this case: h^2 = H^2 = {H_sq_1}. An ideal PGS would explain {max_R_sq_PGS_1*100}% of variance.")
print("  A. PGS can't explain > 50%. (True)")
print("  B. PGS will approach 50%. (True)")
print("  C. PGS will *not* approach 50%. (False)")

# Scenario 2: Non-additive genetic variance exists.
print("\nSCENARIO 2: Non-additive genetic effects exist.")
h_sq_2 = 0.3  # Example value for additive variance component
H_sq_2 = 0.5  # Total genetic variance is still 0.5
max_R_sq_PGS_2 = h_sq_2
print(f"In this case: h^2 = {h_sq_2} < H^2 = {H_sq_2}. An ideal PGS would explain {max_R_sq_PGS_2*100}% of variance.")
print("  A. PGS can't explain > 50%. (True)")
print("  B. PGS will approach 50%. (False)")
print("  C. PGS will *not* approach 50%. (True)")

# Step 4: Conclusion
print("\n--- Final Conclusion ---")
print("A statement is 'necessarily true' only if it holds in ALL possible scenarios.")
print(" - Statement A is TRUE in both scenarios.")
print(" - Statement B is FALSE in Scenario 2.")
print(" - Statement C is FALSE in Scenario 1.")
print("\nTherefore, statement A is the only one that is necessarily true.")
print("The conclusion comes from the definitional equation:")
# This final print statement fulfills the requirement to output the final equation.
# Note: Redirecting output back to normal stdout for the final line.
_print(f"The variance explained by the polygenic score (R^2_PGS) must be less than or equal to the broad-sense heritability (H^2), so: R^2_PGS <= {H_sq}")

<<<A>>>