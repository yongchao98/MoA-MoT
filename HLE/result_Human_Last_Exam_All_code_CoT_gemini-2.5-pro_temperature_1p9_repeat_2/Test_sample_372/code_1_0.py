import math

# Part 1: Calculate the initial percentage of trimethylated sites
# --- Given values ---
final_prop_me3 = 0.1104  # 11.04%
decay_rate_per_hour = 0.10 # 10% or 0.1
time_hours = 10          # 10 hours

# --- Calculation ---
# The formula is P_final = P_initial * (1 - rate)^t
# We solve for P_initial = P_final / (1 - rate)^t
initial_prop_me3 = final_prop_me3 / ((1 - decay_rate_per_hour)**time_hours)
initial_percentage_me3 = initial_prop_me3 * 100

print("--- Part 1: Initial Trimethylation Percentage ---")
print("The calculation for the initial percentage of trimethylated sites is:")
print(f"Initial Percentage = (Final Percentage) / (1 - Decay Rate)^Time")
print(f"Initial Percentage = {final_prop_me3 * 100}% / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"The initial percentage of sites that were trimethylated is: {initial_percentage_me3:.2f}%")
print("\n" + "="*50 + "\n")


# Part 2: Determine the impact on gene expression
# --- Given values ---
initial_expression_rpkm = 200.0
methylation_decrease = 0.10  # 10% decrease

# --- Calculation ---
# Due to the linear relationship, a 10% decrease in methylation proportion
# leads to a 10% decrease in gene expression.
# Final Expression = Initial Expression * (1 - Decrease)
final_expression_rpkm = initial_expression_rpkm * (1 - methylation_decrease)

print("--- Part 2: Impact on Gene Expression ---")
print("The calculation for the new gene expression level is:")
print(f"Final Expression = Initial Expression * (1 - Methylation Decrease)")
print(f"Final Expression = {initial_expression_rpkm} RPKM * (1 - {methylation_decrease})")
print(f"The new gene expression level is: {final_expression_rpkm:.2f} RPKM")
