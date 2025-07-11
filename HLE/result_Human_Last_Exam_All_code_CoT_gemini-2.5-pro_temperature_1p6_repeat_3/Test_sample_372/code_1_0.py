import math

# --- Part 1: Calculate the initial percentage of trimethylated sites ---

# Given values for the first part of the problem
turnover_rate = 0.10  # 10% per hour
time_h = 10  # 10 hours
final_proportion_decimal = 0.1104  # 11.04%

# The decay process is modeled by: Final = Initial * e^(-rate * time)
# To find the initial proportion, we rearrange the formula to: Initial = Final * e^(rate * time)

# Calculate the initial proportion
initial_proportion_decimal = final_proportion_decimal * math.exp(turnover_rate * time_h)
initial_percentage = initial_proportion_decimal * 100

print("--- Part 1: Finding the Initial Methylation Percentage ---")
print("The formula to find the initial proportion A(0) from the final proportion A(t) is: A(0) = A(t) * e^(k*t)")
print(f"Using the given values:")
print(f"Initial Proportion = {final_proportion_decimal} * e^({turnover_rate} * {time_h})")
print(f"Initial Proportion = {final_proportion_decimal} * e^({turnover_rate * time_h})")
print(f"Initial Proportion = {initial_proportion_decimal:.4f}")
print(f"Initial Percentage = {initial_proportion_decimal * 100:.2f}%")
print(f"\nResult: The percentage of sites trimethylated at the start was {initial_percentage:.0f}%.")
print("-" * 50)


# --- Part 2: Determine the impact on gene expression ---

# Given values for the second part of the problem
initial_expression = 200  # RPKM
proportional_decrease = 0.10  # A 10% decrease

# The initial proportion from Part 1 is ~30%
initial_proportion_for_part2 = round(initial_proportion_decimal) # Use the result from Part 1

# Assuming a linear relationship (Expression = m * Proportion), a relative decrease
# in proportion leads to the same relative decrease in expression.
new_expression_level = initial_expression * (1 - proportional_decrease)

print("\n--- Part 2: Determining the Impact on Gene Expression ---")
print("Assuming a linear relationship, a 10% decrease in methylation proportion will cause a 10% decrease in gene expression.")
print("The new gene expression level is calculated as follows:")
print(f"New Expression = Initial Expression * (1 - Proportional Decrease)")
print(f"New Expression = {initial_expression} * (1 - {proportional_decrease})")
print(f"New Expression = {initial_expression} * {1 - proportional_decrease}")
print(f"New Expression = {new_expression_level:.0f}")
print(f"\nResult: The target gene expression would decrease to {new_expression_level:.0f} RPKM.")

# The final answer is the new gene expression level.
# print(f"\n<<<180>>>")