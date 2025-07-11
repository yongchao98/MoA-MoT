import math

# --- Part 1: Calculate the initial percentage of trimethylated sites ---

print("--- Part 1: Calculating the initial percentage of trimethylated sites ---")

# Given values
final_proportion_pct = 11.04  # in percent
time_hours = 10
decay_rate_per_hour = 0.10  # 10% per hour

# Convert percentage to a decimal for calculation
final_proportion = final_proportion_pct / 100

# The formula for exponential decay is P(t) = P(0) * (1 - r)^t
# We need to solve for the initial proportion, P(0).
# P(0) = P(t) / (1 - r)^t
initial_proportion = final_proportion / ((1 - decay_rate_per_hour) ** time_hours)

# Convert the result back to a percentage
initial_proportion_pct = initial_proportion * 100

print(f"The proportion of H3K4me3 sites after {time_hours} hours is {final_proportion_pct}%.")
print(f"The decay rate is {decay_rate_per_hour*100}% per hour.")
print("\nTo find the initial percentage (P0), we use the formula: P0 = P(t) / (1 - r)^t")
print(f"P0 = {final_proportion_pct}% / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"P0 = {final_proportion_pct}% / ({(1 - decay_rate_per_hour):.1f})^{time_hours}")
print(f"P0 = {final_proportion_pct}% / {(1 - decay_rate_per_hour)**time_hours}")
print(f"The initial percentage of trimethylated sites was: {initial_proportion_pct:.2f}%\n")


# --- Part 2: Determine the impact on gene expression ---

print("--- Part 2: Determining the impact on gene expression ---")

# Given values
initial_expression_rpkm = 200
proportion_decrease_pct = 10 # 10% total decrease

# Convert percentage to a decimal
proportion_decrease = proportion_decrease_pct / 100

# Calculate the new expression level
# New Expression = Initial Expression * (1 - Proportion Decrease)
final_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)

print(f"The initial average gene expression is {initial_expression_rpkm} RPKM.")
print(f"The proportion of H3K4me3 sites decreases by {proportion_decrease_pct}%.")
print("\nAssuming a linear relationship, the gene expression level will also decrease by the same proportion.")
print(f"New Expression = Initial Expression * (1 - {proportion_decrease})")
print(f"New Expression = {initial_expression_rpkm} RPKM * {1 - proportion_decrease}")
print(f"The new gene expression level is: {final_expression_rpkm} RPKM")