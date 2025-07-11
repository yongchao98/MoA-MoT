import math

# Part 1: Calculate the initial percentage of H3K4me3 sites

# Define the known values from the problem description
final_percentage_p1 = 11.04  # The percentage of H3K4me3 sites after 10 hours
decay_rate_per_hour = 0.10   # The rate of turnover from H3K4me3 to H3K4me2 per hour
time_hours = 10              # The duration of the period in hours

# The process follows an exponential decay model: Final = Initial * (1 - rate)^time
# We need to find the Initial percentage (P0).
# The formula is rearranged to: Initial = Final / (1 - rate)^time

# Perform the calculation for the initial percentage
decay_factor = (1 - decay_rate_per_hour) ** time_hours
initial_percentage = final_percentage_p1 / decay_factor

# Print the step-by-step calculation for Part 1
print("Part 1: Initial Percentage of H3K4me3 Sites")
print("-------------------------------------------------")
print("The calculation for the initial percentage (P0) is based on the exponential decay formula:")
print(f"P0 = Final_Percentage / (1 - Decay_Rate)^Time")
print(f"P0 = {final_percentage_p1} / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"P0 = {final_percentage_p1} / {decay_factor}")
print(f"The initial percentage of trimethylated sites was: {initial_percentage:.2f}%\n")


# Part 2: Determine the impact on gene expression

# Define the known values for the second part of the problem
initial_expression_rpkm = 200  # Average expression level in RPKM
total_decrease_in_proportion = 0.10 # The total decrease in H3K4me3 proportion over 10 hours

# The relationship between H3K4me3 proportion and gene expression is linear.
# Therefore, a 10% decrease in the methylation proportion leads to a 10% decrease in gene expression.

# Calculate the final gene expression level
final_expression = initial_expression_rpkm * (1 - total_decrease_in_proportion)
change_in_expression = initial_expression_rpkm - final_expression

# Print the step-by-step calculation for Part 2
print("Part 2: Impact on Gene Expression")
print("-------------------------------------------------")
print("Given a linear relationship, the new expression level is calculated as follows:")
print(f"Final_Expression = Initial_Expression * (1 - Total_Decrease_In_Proportion)")
print(f"Final_Expression = {initial_expression_rpkm} * (1 - {total_decrease_in_proportion})")
print(f"Final_Expression = {initial_expression_rpkm} * {1 - total_decrease_in_proportion}")
print(f"The final gene expression level is: {final_expression:.2f} RPKM.")
print(f"This represents a decrease of {change_in_expression:.2f} RPKM from the initial level.")
