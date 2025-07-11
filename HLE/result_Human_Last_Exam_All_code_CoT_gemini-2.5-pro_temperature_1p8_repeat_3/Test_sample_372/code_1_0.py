import math

# Part 1: Calculate the initial percentage of H3K4me3 sites

# Given values
p_final_percent = 11.04  # Final percentage after 10 hours
p_final = p_final_percent / 100.0
decay_rate_per_hour = 0.10  # 10% per hour turnover rate
time_hours = 10

# The formula for exponential decay is P(t) = P0 * (1 - k)^t
# We rearrange to solve for the initial percentage, P0 = P(t) / (1 - k)^t

# Calculation
decay_base = 1 - decay_rate_per_hour
decay_factor = math.pow(decay_base, time_hours)
p_initial = p_final / decay_factor
p_initial_percent = p_initial * 100

print("Part 1: Initial percentage of trimethylated sites")
print("-" * 55)
print("The calculation is based on the exponential decay model.")
print("The equation is: Final Percentage = Initial Percentage * (1 - Decay Rate)^Time")
print("We solve for the Initial Percentage (P0):")
print(f"P0 = {p_final_percent}% / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"P0 = {p_final_percent}% / ({decay_base})^{time_hours}")
print(f"P0 = {p_final_percent}% / {decay_factor:.4f}")
print(f"The calculated initial percentage of H3K4me3 sites is: {p_initial_percent:.2f}%")
print("\n" + "="*55 + "\n")

# Part 2: Determine the impact on gene expression

# Given values
initial_expression_rpkm = 200  # Initial expression level in RPKM
proportion_decrease = 0.10  # 10% decrease in H3K4me3 sites

# Based on the assumption of a linear relationship, a 10% decrease in the
# cause (methylation) leads to a 10% decrease in the effect (gene expression).

# Calculation
final_expression_factor = 1 - proportion_decrease
final_expression_rpkm = initial_expression_rpkm * final_expression_factor

print("Part 2: Impact on target gene expression")
print("-" * 55)
print("This calculation assumes a direct linear relationship between H3K4me3 proportion and gene expression.")
print("A 10% decrease in methylation proportion will cause a 10% decrease in gene expression.")
print("The equation is: New Expression = Initial Expression * (1 - Proportion Decrease)")
print(f"New Expression = {initial_expression_rpkm} RPKM * (1 - {proportion_decrease})")
print(f"New Expression = {initial_expression_rpkm} RPKM * {final_expression_factor}")
print(f"The resulting gene expression level is: {final_expression_rpkm:.2f} RPKM.")

<<<180.0>>>