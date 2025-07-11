import math

# Part 1: Calculate the initial percentage of H3K4me3 sites

# Define the given variables
final_percentage_p_t = 11.04  # The percentage of H3K4me3 sites after 10 hours
turnover_rate_k = 0.10      # The rate of turnover from H3K4me3 to H3K4me2 per hour
time_t = 10                 # The total time in hours

# The formula for exponential decay is P(t) = P0 * (1 - k)^t
# We need to solve for the initial percentage, P0.
# Rearranging the formula gives: P0 = P(t) / (1 - k)^t

# Calculate the initial percentage
decay_factor = (1 - turnover_rate_k) ** time_t
initial_percentage_p0 = final_percentage_p_t / decay_factor

print("--- Part 1: Initial Percentage Calculation ---")
print("To find the initial percentage of trimethylated sites (P0), we use the exponential decay formula.")
print(f"The equation is: P0 = P(t) / (1 - k)^t")
print(f"Plugging in the values: P0 = {final_percentage_p_t} / (1 - {turnover_rate_k})^{time_t}")
print(f"The calculation gives: P0 = {initial_percentage_p0:.2f}%")
print(f"Therefore, the percentage of sites trimethylated at the start was {initial_percentage_p0:.2f}%.")
print("-" * 45)

# Part 2: Calculate the impact on gene expression

# Define the given variables
initial_expression_level = 200  # Average expression in RPKM
methylation_decrease = 0.10     # 10% decrease in the proportion of H3K4me3 sites

# Assuming a direct linear relationship, the expression level will also decrease by the same percentage.
new_expression_level = initial_expression_level * (1 - methylation_decrease)
expression_reduction = initial_expression_level - new_expression_level

print("\n--- Part 2: Gene Expression Impact Calculation ---")
print("Assuming a linear relationship, a 10% decrease in H3K4me3 sites causes a 10% decrease in gene expression.")
print(f"Initial Expression Level: {initial_expression_level} RPKM")
print("To find the new expression level, we reduce the initial level by 10%.")
print(f"New Expression = {initial_expression_level} * (1 - {methylation_decrease})")
print(f"The new expression level is {new_expression_level:.2f} RPKM.")
print(f"\nThe impact is a reduction of {expression_reduction:.2f} RPKM, leading to a new level of {new_expression_level:.2f} RPKM.")
print("-" * 45)
