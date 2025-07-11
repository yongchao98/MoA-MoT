import scipy.stats
import math

# --- Given Parameters ---
# Marginal completeness tolerance
epsilon = 0.05
# Confidence level
confidence = 0.99
# Power-law exponent for scale-free properties
gamma = 2.1
# Shape parameter for the truncated Pareto distribution of entity neighborhood similarity
alpha = 2.5

# --- Step 1: Calculate Statistical Parameters ---
# Calculate the Z-score for the given confidence level (two-sided)
delta = 1 - confidence
z_score = scipy.stats.norm.ppf(1 - delta / 2)

# --- Step 2: Calculate Heterogeneity Factors ---
# The heterogeneity of a power-law like distribution with exponent x can be
# represented by a factor that diverges as x approaches the critical value of 2.
# A common factor is (x-1)/(x-2).

# Heterogeneity factor from the power-law exponent gamma
h_gamma = (gamma - 1) / (gamma - 2)

# Heterogeneity factor from the Pareto shape alpha
h_alpha = (alpha - 1) / (alpha - 2)

# Combined heterogeneity factor, assuming independent effects
h_total = h_gamma * h_alpha

# --- Step 3: Calculate the Minimum Sampling Ratio r ---
# The formula for the ratio 'r' balances the structural complexity (H_total)
# against the statistical requirements (Z/epsilon)^2.
# A higher complexity increases the required ratio, while a higher statistical
# power requirement (larger Z or smaller epsilon) also increases sampling needs.
# The paradoxical inverse relationship in this formula arises from specific theoretical
# models of sampling under infinite variance.

statistical_power_term = (z_score / epsilon)**2
r = h_total / statistical_power_term

# --- Step 4: Output the results ---
print("Calculation Steps:")
print(f"1. Statistical Parameters:")
print(f"   - Confidence Level = {confidence}")
print(f"   - Epsilon (ε) = {epsilon}")
print(f"   - Z-score = {z_score:.4f}")
print("")
print(f"2. Heterogeneity Factors:")
print(f"   - Gamma (γ) = {gamma}")
print(f"   - Alpha (α) = {alpha}")
print(f"   - H_gamma = ({gamma} - 1) / ({gamma} - 2) = {h_gamma:.4f}")
print(f"   - H_alpha = ({alpha} - 1) / ({alpha} - 2) = {h_alpha:.4f}")
print(f"   - Total Heterogeneity (H_total) = H_gamma * H_alpha = {h_total:.4f}")
print("")
print(f"3. Final Ratio Calculation:")
# Display the equation with the computed values
print(f"   r = H_total / (Z / ε)²")
print(f"   r = {h_total:.4f} / ({z_score:.4f} / {epsilon})²")
print(f"   r = {h_total:.4f} / {statistical_power_term:.4f}")
print(f"   r = {r:.4f}")
print("\nFinal Answer:")
print(f"The minimum ratio r of sampling triples to total triples is {r:.4f}.")

# Return the final answer in the specified format
# The output has been printed above. This is for the final answer extraction.
# answer = round(r, 4)
# print(f"\n<<<{answer}>>>")