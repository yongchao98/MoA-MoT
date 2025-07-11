import math

# --- Part 1: Calculate the initial percentage of trimethylated sites ---

# Given values
final_proportion_p1 = 11.04  # in percent
decay_rate = 0.10          # 10% per hour
time = 10                  # in hours

# The formula for exponential decay is P(t) = P0 * (1 - k)^t
# We need to solve for the initial proportion, P0.
# P0 = P(t) / (1 - k)^t
base = 1 - decay_rate
initial_proportion = final_proportion_p1 / (base ** time)

print("--- Part 1: Initial Percentage Calculation ---")
print("The problem describes an exponential decay of H3K4me3 sites.")
print(f"To find the initial percentage, we use the formula: Initial = Final / (1 - Rate)^Time")
print(f"Initial Percentage = {final_proportion_p1} / (1 - {decay_rate})^{time}")
print(f"Initial Percentage = {final_proportion_p1} / ({base:.1f})^{time}")
print(f"Initial Percentage = {final_proportion_p1} / {base**time:.4f}")
print(f"\nThe initial percentage of H3K4me3 sites was: {initial_proportion:.2f}%\n")


# --- Part 2: Determine the impact on gene expression ---

# Given values and assumptions
max_expression = 200  # RPKM for 100% methylation
max_proportion = 100  # 100%
relative_decrease = 0.10 # 10%

# Establish the linear relationship: E = m * P
# Assuming E=0 at P=0, and E=200 at P=100.
# m = 200 / 100 = 2
slope = max_expression / max_proportion

# Calculate the initial expression level using the result from Part 1
initial_expression = slope * initial_proportion

# Calculate the new proportion after a 10% relative decrease
final_proportion_p2 = initial_proportion * (1 - relative_decrease)

# Calculate the new expression level
final_expression = slope * final_proportion_p2

# Calculate the change in expression
expression_change = final_expression - initial_expression

print("--- Part 2: Gene Expression Impact Analysis ---")
print("We assume a linear relationship between methylation and expression: E = m * P.")
print(f"With {max_expression} RPKM at {max_proportion}% methylation, the slope (m) is {slope:.1f}.")
print(f"\nInitial State:")
print(f"The initial methylation of {initial_proportion:.2f}% corresponds to an expression level of:")
print(f"Expression_initial = {slope:.1f} * {initial_proportion:.2f} = {initial_expression:.2f} RPKM")

print(f"\nFinal State (after a {relative_decrease*100}% decrease in methylation):")
print(f"New Methylation = {initial_proportion:.2f}% * (1 - {relative_decrease}) = {final_proportion_p2:.2f}%")
print(f"This new methylation level corresponds to an expression level of:")
print(f"Expression_final = {slope:.1f} * {final_proportion_p2:.2f} = {final_expression:.2f} RPKM")

print(f"\nImpact on Gene Expression:")
print(f"Change = Final Expression - Initial Expression")
print(f"Change = {final_expression:.2f} - {initial_expression:.2f} = {expression_change:.2f} RPKM")
print(f"The impact is a decrease of {-expression_change:.2f} RPKM in the target gene expression.")

# Final answer in the required format
final_answer = -expression_change
# <<<6.33>>>