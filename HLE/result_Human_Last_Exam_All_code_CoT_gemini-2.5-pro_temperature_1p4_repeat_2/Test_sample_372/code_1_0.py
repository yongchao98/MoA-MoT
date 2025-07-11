import math

# Part 1: Calculate the initial percentage of trimethylated sites.
# Define the known values from the problem.
p_t = 11.04  # Percentage of H3K4me3 after 10 hours.
k = 0.10     # Turnover rate per hour (10%).
t = 10       # Time in hours.

# The process is modeled by the exponential decay formula: P(t) = P(0) * e^(-k*t).
# We rearrange the formula to solve for the initial percentage, P(0).
# P(0) = P(t) / e^(-k*t)
p_0 = p_t / math.exp(-k * t)

print("Part 1: Calculating the initial percentage of trimethylated sites")
print("----------------------------------------------------------------")
print(f"The model for the decay of H3K4me3 is: P(t) = P(0) * e^(-k*t)")
print(f"Given values: P({t}) = {p_t}%, rate k = {k}, time t = {t} hours.")
print("To find the initial percentage P(0), we rearrange the formula:")
print(f"P(0) = P({t}) / e^(-k*t)")
print("Plugging in the numbers:")
print(f"P(0) = {p_t} / e^(-{k} * {t})")
print(f"P(0) = {round(p_0, 2)}%")
print(f"Therefore, the initial percentage of trimethylated sites was {round(p_0, 2)}%.\n")


# Part 2: Determine the impact on gene expression.
# This is a hypothetical scenario based on the result from Part 1.
initial_expression = 200  # Initial gene expression in RPKM.
initial_proportion = p_0    # The calculated initial methylation percentage.
decrease_percentage = 10    # The hypothetical 10% decrease.

# With a linear relationship, a 10% relative decrease in methylation
# leads to a 10% relative decrease in gene expression.
new_expression = initial_expression * (1 - (decrease_percentage / 100.0))

print("Part 2: Calculating the impact on target gene expression")
print("---------------------------------------------------------")
print(f"This is a hypothetical scenario assuming a {decrease_percentage}% decrease in H3K4me3 from the initial state.")
print(f"Initial H3K4me3 proportion: {round(initial_proportion, 2)}%")
print(f"Initial gene expression: {initial_expression} RPKM")
print("Assuming a linear relationship, a 10% decrease in methylation causes a 10% decrease in expression.")
print("The calculation for the new expression level is:")
print(f"New Expression = Initial Expression * (1 - ({decrease_percentage} / 100))")
print(f"New Expression = {initial_expression} * (1 - {decrease_percentage / 100.0})")
print(f"New Expression = {int(new_expression)} RPKM")
print(f"\nThe impact on the target gene expression is a new level of {int(new_expression)} RPKM.")

<<<180>>>