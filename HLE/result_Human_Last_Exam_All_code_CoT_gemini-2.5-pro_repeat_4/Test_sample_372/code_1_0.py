import math

# --- Given parameters ---
# Part 1 parameters
final_percentage = 11.04  # P(t) in percent
decay_rate = 0.10         # k, 10% per hour
time = 10                 # t, in hours

# Part 2 parameters
initial_expression = 200  # RPKM
proportion_decrease = 0.10 # 10%

# --- Part 1: Calculate the initial percentage of H3K4me3 ---
print("--- Part 1: Calculating the Initial H3K4me3 Percentage ---")
print("The decay follows the formula: P(t) = P0 * (1 - k)^t")
print("We need to solve for the initial percentage (P0).")
print(f"P0 = P(t) / (1 - k)^t")

# Calculate the decay factor
decay_factor = (1 - decay_rate) ** time

# Calculate the initial percentage
initial_percentage = final_percentage / decay_factor

print(f"P0 = {final_percentage} / (1 - {decay_rate})^{time}")
print(f"P0 = {final_percentage} / {decay_factor:.6f}")
print(f"The initial percentage of H3K4me3 sites is {initial_percentage:.2f}%.")
print("-" * 55)

# --- Part 2: Determine the impact on gene expression ---
print("\n--- Part 2: Calculating the Impact on Gene Expression ---")
print("Assuming a direct linear relationship between H3K4me3 proportion and gene expression.")
print("A 10% decrease in methylation proportion will cause a 10% decrease in expression.")

# Calculate the decrease in expression
expression_decrease = initial_expression * proportion_decrease

# Calculate the final expression
final_expression = initial_expression - expression_decrease

print(f"Initial Gene Expression = {initial_expression} RPKM")
print(f"Decrease in Expression = {initial_expression} * {proportion_decrease}")
print(f"Final Gene Expression = {initial_expression} - ({initial_expression} * {proportion_decrease})")
print(f"Final Gene Expression = {initial_expression} - {expression_decrease}")
print(f"The new gene expression level is {final_expression} RPKM.")
print("-" * 55)
