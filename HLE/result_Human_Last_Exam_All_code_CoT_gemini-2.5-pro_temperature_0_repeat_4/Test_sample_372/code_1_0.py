import math

# Part 1: Calculate the initial percentage of trimethylated sites.

# Given values from the problem
p_final_percent = 11.04  # Percentage of H3K4me3 sites after 10 hours
decay_rate_per_hour = 0.10 # 10% per hour rate of turnover
time_hours = 10

# The process follows an exponential decay model: P(t) = P0 * (1 - r)^t
# We need to find the initial percentage, P0, by rearranging the formula: P0 = P(t) / (1 - r)^t

# Perform the calculation for P0
decay_factor = 1 - decay_rate_per_hour
total_decay_multiplier = math.pow(decay_factor, time_hours)
p_initial_percent = p_final_percent / total_decay_multiplier

print("--- Part 1: Initial Methylation Percentage Calculation ---")
print("The formula for exponential decay is: P(t) = P0 * (1 - r)^t")
print(f"To find the initial percentage (P0), we use: P0 = P(t) / (1 - r)^t")
print(f"P0 = {p_final_percent} / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"P0 = {p_final_percent} / {total_decay_multiplier:.6f}")
print(f"Result: The initial percentage of trimethylated sites (P0) was {p_initial_percent:.2f}%\n")


# Part 2: Determine the impact on gene expression.

# Given values from the problem
initial_expression_rpkm = 200 # Average expression level in RPKM
proportion_decrease = 0.10 # A 10% decrease in H3K4me3 sites

# A linear relationship means that gene expression is directly proportional to the methylation proportion.
# A 10% decrease in proportion means the new proportion is 90% of the original.
# Therefore, the new expression level will be 90% of the original expression level.

# Perform the calculation for the new expression level
final_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)

print("--- Part 2: Gene Expression Impact Calculation ---")
print("Assuming a linear relationship, a 10% decrease in the methylation proportion results in a 10% decrease in gene expression.")
print("The calculation is: New Expression = Initial Expression * (1 - Percentage Decrease)")
print(f"New Expression = {initial_expression_rpkm} * (1 - {proportion_decrease})")
print(f"New Expression = {initial_expression_rpkm} * {1 - proportion_decrease}")
print(f"Result: The impacted gene expression level is {final_expression_rpkm:.2f} RPKM")
