import math

# Part 1: Calculate the initial percentage of trimethylated sites

# Known values from the problem statement
final_proportion_perc = 11.04
decay_rate_per_hour = 0.10
time_in_hours = 10

# Convert percentage to a decimal for calculation
final_proportion = final_proportion_perc / 100.0

# Formula: P_initial = P_final / (1 - rate)^time
initial_proportion = final_proportion / math.pow((1 - decay_rate_per_hour), time_in_hours)

# Convert the result back to a percentage
initial_proportion_perc = initial_proportion * 100

print("Part 1: Calculation of the initial H3K4me3 percentage")
print("======================================================")
print(f"The final percentage of H3K4me3 after {time_in_hours} hours is {final_proportion_perc}%.")
print(f"The decay rate is {decay_rate_per_hour*100}% per hour.")
print(f"To find the initial percentage (P_initial), we use the equation: P_initial = P_final / (1 - rate)^time")
print(f"Calculation: P_initial = {final_proportion_perc}% / (1 - {decay_rate_per_hour})^{time_in_hours}")
print(f"The initial percentage of H3K4me3 sites was {initial_proportion_perc:.2f}%.\n")


# Part 2: Determine the impact on gene expression

# Known values from the problem statement
initial_expression_rpkm = 200
proportion_decrease_perc = 10.0

# Due to the linear relationship, a 10% decrease in the methylation proportion
# results in a 10% decrease in gene expression.
decrease_factor = proportion_decrease_perc / 100.0
final_expression_rpkm = initial_expression_rpkm * (1 - decrease_factor)

print("Part 2: Calculation of the new gene expression level")
print("=====================================================")
print(f"The initial gene expression is {initial_expression_rpkm} RPKM.")
print(f"Assuming a {proportion_decrease_perc}% decrease in the H3K4me3 proportion.")
print(f"To find the new expression level (E_final), we use the equation: E_final = E_initial * (1 - {proportion_decrease_perc}%)")
print(f"Calculation: E_final = {initial_expression_rpkm} RPKM * (1 - {decrease_factor})")
print(f"The new target gene expression level would be {final_expression_rpkm:.2f} RPKM.")

<<<180.00>>>