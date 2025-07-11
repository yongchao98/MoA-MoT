import math

# Part 1: Calculate the initial percentage of H3K4me3 sites

# Define given values for the decay process
p_final_percent = 11.04  # Final proportion in percent after 10 hours
k = 0.10  # Decay rate per hour
t = 10  # Time in hours

print("Part 1: Calculating the initial percentage of H3K4me3.")
print("The process is modeled by exponential decay: P(t) = P(0) * e^(-k*t)")
print(f"We are given P({t}) = {p_final_percent}%, rate k = {k}/hour, and time t = {t} hours.")
print("We need to solve for the initial percentage, P(0).")
print(f"The equation is: {p_final_percent} = P(0) * e^(-{k} * {t})")

# To solve for P(0), we rearrange the formula: P(0) = P(t) * e^(k*t)
p_initial_calculated = p_final_percent * math.exp(k * t)

# The result is very close to 30, so we round it.
p_initial_rounded = round(p_initial_calculated)

print("\nSolving for P(0):")
print(f"P(0) = {p_final_percent} * e^({k} * {t})")
print(f"P(0) = {p_final_percent} * e^{k*t}")
print(f"P(0) = {p_initial_calculated:.2f}%")
print(f"Based on this result, we can assume the initial percentage was {p_initial_rounded}%.")

print("\n" + "-"*50 + "\n")

# Part 2: Calculate the new gene expression level

# Define given values for the gene expression calculation
e_initial_rpkm = 200  # Initial expression level in RPKM
p_initial_percent = p_initial_rounded # Use the rounded initial proportion
proportion_decrease = 0.10  # The proportion of H3K4me3 decreases by 10% (0.10)

print("Part 2: Calculating the impact on target gene expression.")
print(f"The initial state is {p_initial_percent}% H3K4me3 with an expression of {e_initial_rpkm} RPKM.")
print(f"The H3K4me3 proportion then decreases by {proportion_decrease*100}%.")

# Calculate the new proportion of H3K4me3
decrease_amount = p_initial_percent * proportion_decrease
p_new_percent = p_initial_percent - decrease_amount

print("\nCalculating the new H3K4me3 proportion:")
print(f"New Proportion = Initial Proportion * (1 - Decrease Fraction)")
print(f"New Proportion = {p_initial_percent}% * (1 - {proportion_decrease})")
print(f"New Proportion = {p_new_percent}%")


# Assuming a linear relationship E = m * P, where E is expression and P is proportion.
# First, find the proportionality constant m from the initial state.
# m = E_initial / P_initial
# Then calculate the new expression E_new = m * P_new
e_new_rpkm = (e_initial_rpkm / p_initial_percent) * p_new_percent

print("\nCalculating the new gene expression level:")
print("Based on the linear relationship: New Expression = (Initial Expression / Initial Proportion) * New Proportion")
print(f"New Expression = ({e_initial_rpkm} RPKM / {p_initial_percent}%) * {p_new_percent}%")
print(f"New Expression = {e_new_rpkm:.0f} RPKM")

print("\nThe impact on gene expression is a decrease to the new level.")
print(f"Final calculated gene expression: {e_new_rpkm:.0f} RPKM")

<<<180>>>