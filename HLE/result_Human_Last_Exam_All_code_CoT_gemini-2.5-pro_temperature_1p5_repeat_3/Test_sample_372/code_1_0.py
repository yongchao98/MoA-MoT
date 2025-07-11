import math

# --- Part 1: Calculate the initial percentage of H3K4me3 sites ---
# The decay of methylation is modeled using the continuous exponential decay formula:
# P(t) = P(0) * e^(-k*t)
# where P(t) is the proportion at time t, P(0) is the initial proportion,
# k is the rate constant, and t is time.

# Given values for Part 1
p_final = 11.04    # Final proportion in percent
rate = 0.10        # 10% per hour decay rate
time = 10          # 10 hours

# Solving for the initial proportion P_initial
p_initial = p_final / math.exp(-rate * time)

print("--- Part 1: Initial Methylation Calculation ---")
print("To find the initial percentage of trimethylated sites (P_initial), we use the formula:")
print("P_initial = P_final / e^(-rate * time)")
print(f"Plugging in the values:")
print(f"P_initial = {p_final}% / e^(-{rate} * {time})")
print(f"The initial percentage was {p_initial:.2f}%.\n")


# --- Part 2: Determine the impact on gene expression ---
# We are given a linear relationship between methylation proportion and gene expression.
# We'll assume direct proportionality: Expression = Constant * Proportion.
# If the initial expression is E_initial at proportion P_initial, then E_initial = k * P_initial.
# A 10% decrease in proportion means the new proportion is P_final = P_initial * (1 - 0.10).
# The new expression E_final = k * P_final = k * (P_initial * 0.9) = (k * P_initial) * 0.9 = E_initial * 0.9.
# Thus, a 10% relative decrease in methylation proportion leads to a 10% decrease in gene expression.

# Given values for Part 2
e_initial = 200    # Initial gene expression in RPKM
decrease = 0.10    # 10% decrease

# Calculation for Part 2
e_final = e_initial * (1 - decrease)

print("--- Part 2: Gene Expression Impact Calculation ---")
print("To find the new gene expression level, we calculate the effect of a 10% decrease:")
print(f"Final Expression = Initial Expression * (1 - Decrease Percentage)")
print(f"Final Expression = {e_initial} RPKM * (1.0 - {decrease}) = {e_final} RPKM")

print(f"\n<<<{e_final}>>>")