import math

# Part 1: Calculating the initial percentage of H3K4me3 sites.
# We model the turnover as a continuous exponential decay.
print("--- Part 1: Calculating the initial percentage of H3K4me3 sites ---")

# Given values for Part 1
rate = 0.10  # 10% per hour
time = 10.0  # 10 hours
final_proportion_part1 = 0.1104  # 11.04%

# The formula for exponential decay is P(t) = P(0) * e^(-r*t).
# We rearrange it to solve for the initial proportion P(0).
# P(0) = P(t) / e^(-r*t)
initial_proportion = final_proportion_part1 / math.exp(-rate * time)
initial_percentage = initial_proportion * 100

print("The model for H3K4me3 turnover is: P_final = P_initial * e^(-rate * time)")
print("We solve for the initial proportion (P_initial):")
print(f"P_initial = {final_proportion_part1} / e^(-{rate} * {time})")
print(f"P_initial = {initial_proportion:.4f} or {initial_percentage:.2f}%")
print(f"The percentage of sites trimethylated at the start was {initial_percentage:.2f}%.")
print("-" * 50)

# Part 2: Determining the impact on gene expression.
# The relationship between H3K4me3 proportion and gene expression is linear.
print("\n--- Part 2: Determining the impact on gene expression ---")

# Given values for Part 2
initial_expression = 200.0  # RPKM
decrease_in_proportion = 0.10  # 10% decrease

# A linear relationship means Expression = k * Proportion.
# If the proportion decreases by 10%, the expression also decreases by 10%.
new_expression = initial_expression * (1 - decrease_in_proportion)

print("A 10% decrease in the proportion of H3K4me3 sites results in a 10% decrease in gene expression.")
print("The new expression level is calculated as follows:")
print(f"New Expression = Initial Expression * (1 - Proportional Decrease)")
print(f"New Expression = {initial_expression} * (1 - {decrease_in_proportion}) = {new_expression:.1f} RPKM")
print(f"The target gene expression level would be reduced to {new_expression:.1f} RPKM.")

<<<180.0>>>