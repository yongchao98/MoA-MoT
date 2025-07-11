import math

# Part 1: Calculate the initial percentage of H3K4me3

# Given values
p_final_percent = 11.04  # Percentage of H3K4me3 after 10 hours
rate_k = 0.10            # 10% per hour turnover rate
time_t = 10.0            # 10 hours

# The formula for exponential decay is P(t) = P0 * e^(-k*t)
# We rearrange to solve for the initial percentage, P0: P0 = P(t) / e^(-k*t)
p_initial_percent = p_final_percent / math.exp(-rate_k * time_t)
# Round to a reasonable precision
p_initial_percent = round(p_initial_percent, 2)

print("--- Part 1: Initial H3K4me3 Percentage ---")
print(f"The final percentage of H3K4me3 after {time_t} hours is {p_final_percent}%.")
print(f"The turnover rate is {rate_k*100}% per hour.")
print("Using the exponential decay formula: P(t) = P0 * e^(-k*t)")
print(f"We solve for the initial percentage P0: P0 = P(t) / e^(-k*t)")
print(f"P0 = {p_final_percent} / e^(-{rate_k} * {time_t})")
print(f"The initial percentage of H3K4me3 sites (P0) was {p_initial_percent}%.")
print("\n" + "="*50 + "\n")


# Part 2: Determine the impact on gene expression

# Given values
e_initial_rpkm = 200.0   # Initial gene expression level in RPKM
percent_decrease = 0.10  # 10% decrease in the proportion of H3K4me3 sites

# Calculate the new proportion of H3K4me3
p_new_percent = p_initial_percent * (1 - percent_decrease)

# Assuming a linear relationship (direct proportionality):
# New Expression / Initial Expression = New Proportion / Initial Proportion
# New Expression = Initial Expression * (New Proportion / Initial Proportion)
e_new_rpkm = e_initial_rpkm * (p_new_percent / p_initial_percent)

print("--- Part 2: Impact on Gene Expression ---")
print("Assuming a linear relationship between H3K4me3 proportion and gene expression.")
print(f"Initial state: A proportion of {p_initial_percent}% H3K4me3 corresponds to {e_initial_rpkm} RPKM.")
print(f"The proportion of H3K4me3 sites decreases by {percent_decrease*100}%.")
print(f"The new H3K4me3 proportion is {p_initial_percent}% * (1 - {percent_decrease}) = {p_new_percent}%.")
print("Calculating the new gene expression level:")
print(f"New Expression = {e_initial_rpkm} RPKM * ({p_new_percent}% / {p_initial_percent}%)")
print(f"The new gene expression level is {e_new_rpkm} RPKM.")

# The final answer required by the prompt
# This is the new gene expression level, which quantifies the impact
final_answer = e_new_rpkm
# <<<180.0>>>