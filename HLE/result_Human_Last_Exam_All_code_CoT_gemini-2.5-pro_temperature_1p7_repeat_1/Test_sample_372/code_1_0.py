import math

# Part 1: Calculate the initial percentage of trimethylated sites

# Given values
p_final_percent = 11.04  # Percentage of H3K4me3 sites after 10 hours
decay_rate_per_hour = 0.10 # 10% turnover rate per hour
time_hours = 10          # The period is 10 hours

# The formula for exponential decay is P(t) = P0 * (1 - r)^t
# We need to solve for the initial percentage, P0.
# P0 = P(t) / (1 - r)^t
initial_factor = (1 - decay_rate_per_hour)**time_hours
initial_percentage = p_final_percent / initial_factor

print("--- Part 1: Initial Methylation Percentage ---")
print("To find the initial percentage of trimethylated sites (P0), we use the exponential decay formula.")
print(f"The equation is: P0 = P(t) / (1 - r)^t")
print(f"Using the given values:")
print(f"P0 = {p_final_percent}% / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"P0 = {p_final_percent}% / {initial_factor:.4f}")
print(f"The initial percentage of trimethylated sites was {initial_percentage:.2f}%.")
print("-" * 50)


# Part 2: Determine the impact on gene expression

# Given values
initial_expression_rpkm = 200  # RPKM at the initial methylation state
proportion_decrease = 0.10   # 10% decrease in the proportion of H3K4me3 sites

# Assuming a linear relationship, a 10% decrease in methylation proportion
# leads to a 10% decrease in gene expression.
final_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)

print("\n--- Part 2: Impact on Gene Expression ---")
print("Assuming a linear relationship between methylation and gene expression, a 10% decrease in the H3K4me3 proportion will cause a 10% decrease in expression.")
print(f"The calculation for the new expression level (E_final) is:")
print(f"E_final = Initial Expression * (1 - Proportion Decrease)")
print(f"E_final = {initial_expression_rpkm} RPKM * (1 - {proportion_decrease})")
print(f"The target gene expression level would decrease to {final_expression_rpkm:.2f} RPKM.")

# The final answer required is the new expression level.
final_answer = final_expression_rpkm

# This line is for the final answer format and will be hidden in the assistant's output logic
# print(f'<<<{final_answer}>>>')