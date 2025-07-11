import math

# Part 1: Calculate the initial percentage of trimethylated sites.
# Given values
p_final_percentage = 11.04  # Final percentage of H3K4me3 sites
decay_rate = 0.10           # 10% decay rate per hour
time_hours = 10             # 10 hours duration

# The formula for exponential decay is P_final = P_initial * (1 - r)^t
# We solve for P_initial: P_initial = P_final / (1 - r)^t
p_initial_percentage = p_final_percentage / ((1 - decay_rate)**time_hours)

print("--- Part 1: Initial Percentage of H3K4me3 Sites ---")
print("The calculation is based on the exponential decay model.")
print(f"The equation to solve is: {p_final_percentage} = P0 * (1 - {decay_rate})^{time_hours}")
print(f"Rearranging to solve for the initial percentage (P0): P0 = {p_final_percentage} / ({1 - decay_rate})^{time_hours}")
print(f"The initial percentage of sites that were trimethylated is: {p_initial_percentage:.2f}%\n")


# Part 2: Determine the impact on gene expression.
# Given values
initial_expression_rpkm = 200  # Initial gene expression level
methylation_decrease = 0.10    # 10% decrease in H3K4me3 sites

# Due to the linear relationship, a 10% decrease in methylation causes a 10% decrease in expression.
final_expression_rpkm = initial_expression_rpkm * (1 - methylation_decrease)
expression_change_rpkm = initial_expression_rpkm - final_expression_rpkm

print("--- Part 2: Impact on Gene Expression ---")
print("The calculation is based on a linear relationship between methylation and expression.")
print(f"The equation for the new expression level is: New Expression = {initial_expression_rpkm} * (1 - {methylation_decrease})")
print(f"The new gene expression level is: {final_expression_rpkm:.2f} RPKM.")
print(f"This is a decrease of {expression_change_rpkm:.2f} RPKM from the initial level.")

<<<31.66>>>