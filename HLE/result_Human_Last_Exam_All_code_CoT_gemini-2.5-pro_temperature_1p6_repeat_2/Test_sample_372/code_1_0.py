import math

# Part 1: Calculate the initial percentage of H3K4me3

# Given values
P_t = 11.04  # Percentage of H3K4me3 at t=10 hours
rate = 0.10  # 10% decay rate per hour
t = 10       # 10 hours

# The formula for continuous decay is P(t) = P0 * e^(-rate*t)
# We rearrange to solve for the initial percentage, P0
# P0 = P(t) / e^(-rate*t)
P0_calculated = P_t / math.exp(-rate * t)

# The calculated value is very close to 30, so we'll use 30 as the initial percentage.
P0_assumed = round(P0_calculated)

print("--- Part 1: Initial Methylation Percentage ---")
print(f"The calculation for the initial percentage (P0) is based on the formula: P0 = P(t) * e^(k*t)")
print(f"P0 = {P_t}% * e^({rate} * {t})")
print(f"P0 = {P0_calculated:.2f}%")
print(f"Based on this, the initial percentage of trimethylated sites is assumed to be: {P0_assumed}%\n")


# Part 2: Calculate the impact on gene expression

# Given values
GE_initial = 200  # RPKM
P_initial = P0_assumed # Using the result from Part 1
decrease = 0.10   # 10% decrease in proportion

# Calculate the final methylation proportion
P_final = P_initial * (1 - decrease)

# Assuming a linear relationship (direct proportionality), GE_final / P_final = GE_initial / P_initial
# Therefore, GE_final = GE_initial * (P_final / P_initial)
# This simplifies to GE_final = GE_initial * (1 - decrease)
GE_final = GE_initial * (1 - decrease)

print("--- Part 2: Impact on Gene Expression ---")
print(f"With an initial methylation of {P_initial}% corresponding to {GE_initial} RPKM, a {decrease*100}% decrease in methylation leads to a new expression level.")
print("The final gene expression (GE_final) is calculated as: GE_final = GE_initial * (1 - decrease)")
print(f"{GE_final:.1f} RPKM = {GE_initial} RPKM * (1 - {decrease})")
print(f"The new average expression level of the target genes is: {GE_final:.1f} RPKM")

# Final Answer
# Format: <<<Initial percentage: [value]%, Final expression: [value] RPKM>>>
print("\n<<<Initial percentage: 30.0%, Final expression: 180.0 RPKM>>>")