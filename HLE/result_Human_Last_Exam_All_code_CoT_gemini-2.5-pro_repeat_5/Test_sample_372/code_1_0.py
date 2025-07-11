import math

# Part 1: Calculate the initial percentage of H3K4me3 sites

# Given values from the problem
final_proportion_percent = 11.04
decay_rate_per_hour = 0.10
time_hours = 10

# The process is an exponential decay. The formula is:
# Final_Proportion = Initial_Proportion * (1 - Decay_Rate)^Time
# We need to solve for the Initial_Proportion.
# Initial_Proportion = Final_Proportion / (1 - Decay_Rate)^Time

# Perform the calculation
decay_factor = (1 - decay_rate_per_hour) ** time_hours
initial_proportion_percent = final_proportion_percent / decay_factor

print("--- Part 1: Initial H3K4me3 Percentage Calculation ---")
print("To find the initial percentage, we use the exponential decay formula rearranged to solve for the initial amount.")
print(f"Initial Percentage = Final Percentage / (1 - Decay Rate)^Time")
print(f"Initial Percentage = {final_proportion_percent}% / (1 - {decay_rate_per_hour})^{time_hours}")
# Displaying the calculated decay factor before the final division for clarity
print(f"Initial Percentage = {final_proportion_percent}% / {decay_factor:.6f}")
print(f"Result: The initial percentage of trimethylated sites was {initial_proportion_percent:.2f}%.")
print("\n" + "="*60 + "\n")

# Part 2: Determine the impact on target gene expression

# Given values from the problem
initial_expression_rpkm = 200
proportion_decrease = 0.10  # This represents a 10% decrease

# The relationship is linear, so a 10% decrease in H3K4me3 proportion
# leads to a 10% decrease in gene expression.
# New Expression = Initial Expression * (1 - Proportional Decrease)

# Perform the calculation
final_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)

print("--- Part 2: Gene Expression Impact Calculation ---")
print("Assuming a linear relationship, a 10% decrease in the cause (methylation) leads to a 10% decrease in the effect (expression).")
print(f"New Expression Level = Initial Expression * (1 - Proportional Decrease)")
print(f"New Expression Level = {initial_expression_rpkm} RPKM * (1 - {proportion_decrease})")
print(f"New Expression Level = {initial_expression_rpkm} RPKM * {1 - proportion_decrease}")
print(f"Result: The new target gene expression level is {final_expression_rpkm:.2f} RPKM.")
