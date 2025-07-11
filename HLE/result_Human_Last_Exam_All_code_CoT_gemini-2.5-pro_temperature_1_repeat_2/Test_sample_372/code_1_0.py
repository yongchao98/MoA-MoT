import math

# --- Given values ---
# Part 1 values
turnover_rate_per_hour = 0.10  # 10%
time_hours = 10
final_proportion_H3K4me3_percent = 11.04

# Part 2 values
initial_expression_RPKM = 200
proportion_decrease_percent = 10.0

# --- Part 1: Calculate the initial percentage of H3K4me3 sites ---

# Convert percentage to decimal for calculation
final_proportion = final_proportion_H3K4me3_percent / 100

# The decay factor is (1 - turnover_rate_per_hour)^time_hours
decay_factor = (1 - turnover_rate_per_hour) ** time_hours

# Calculate the initial proportion by rearranging the exponential decay formula: P(0) = P(t) / (1-r)^t
initial_proportion = final_proportion / decay_factor

# Convert the result back to a percentage
initial_proportion_percent = initial_proportion * 100

print("--- Part 1: Initial Methylation Calculation ---")
print(f"The calculation for the initial percentage of H3K4me3 sites is:")
print(f"Initial % = Final % / (1 - Rate)^Time")
print(f"Initial % = {final_proportion_H3K4me3_percent:.2f}% / (1 - {turnover_rate_per_hour})^{time_hours}")
print(f"Initial % = {final_proportion_H3K4me3_percent:.2f}% / {decay_factor:.4f}")
print(f"Result: The initial percentage of trimethylated sites was {initial_proportion_percent:.2f}%.")


# --- Part 2: Determine the impact on gene expression ---

# Calculate the new expression level after a 10% decrease in methylation proportion.
# New Expression = Initial Expression * (1 - Percentage Decrease)
decrease_factor = 1 - (proportion_decrease_percent / 100)
new_expression_RPKM = initial_expression_RPKM * decrease_factor

print("\n--- Part 2: Gene Expression Impact Calculation ---")
print(f"The calculation for the new gene expression level is:")
print(f"New Expression = Initial Expression * (1 - Proportion Decrease)")
print(f"New Expression = {initial_expression_RPKM} RPKM * (1 - {proportion_decrease_percent/100})")
print(f"Result: The new gene expression level is {new_expression_RPKM:.0f} RPKM.")
