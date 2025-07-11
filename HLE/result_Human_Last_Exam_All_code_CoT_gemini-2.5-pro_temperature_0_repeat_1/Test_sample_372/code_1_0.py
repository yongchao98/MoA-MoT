import math

# Part 1: Calculate the initial percentage of H3K4me3

# Given values for the first part of the problem
p_final_percentage = 11.04  # Percentage of H3K4me3 after 10 hours
decay_rate_per_hour = 0.10      # 10% decay rate per hour
time_hours = 10        # Time in hours

# The formula for exponential decay is P(t) = P0 * (1 - r)^t
# We rearrange it to solve for the initial percentage P0: P0 = P(t) / (1 - r)^t
denominator = (1 - decay_rate_per_hour) ** time_hours
p_initial_percentage = p_final_percentage / denominator

print("--- Part 1: Initial Percentage Calculation ---")
print("To find the initial percentage of trimethylated sites (P0), we use the exponential decay formula.")
print("\nEquation: P0 = P(t) / (1 - r)^t")
print(f"P(t) = {p_final_percentage}% (final percentage)")
print(f"r = {decay_rate_per_hour} (decay rate)")
print(f"t = {time_hours} hours (time)")
print(f"\nCalculation: P0 = {p_final_percentage} / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"Result: The initial percentage of trimethylated sites was {p_initial_percentage:.2f}%\n")


# Part 2: Calculate the impact on gene expression

# Given values for the second part of the problem
initial_expression_rpkm = 200  # Initial gene expression in RPKM
proportion_decrease = 0.10 # 10% decrease in H3K4me3 sites

# Due to the stated linear relationship, a 10% decrease in methylation proportion
# leads to a 10% decrease in gene expression.
new_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)
expression_change = new_expression_rpkm - initial_expression_rpkm

print("--- Part 2: Gene Expression Impact Calculation ---")
print("Due to the linear relationship, a decrease in methylation proportion causes a proportional decrease in gene expression.")
print("\nEquation: E_new = E_initial * (1 - decrease_fraction)")
print(f"E_initial = {initial_expression_rpkm} RPKM (initial expression)")
print(f"decrease_fraction = {proportion_decrease} (10% decrease)")
print(f"\nCalculation: E_new = {initial_expression_rpkm} * (1 - {proportion_decrease})")
print(f"Result: The new gene expression level is {new_expression_rpkm:.2f} RPKM.")
print(f"This is a decrease of {abs(expression_change):.2f} RPKM from the initial level.")

# The final answer is the new gene expression level.
final_answer = new_expression_rpkm
# <<<180.0>>>