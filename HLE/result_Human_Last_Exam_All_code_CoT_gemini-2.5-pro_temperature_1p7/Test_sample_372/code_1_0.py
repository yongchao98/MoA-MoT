import math

# Part 1: Calculate the initial percentage of trimethylated sites
print("Part 1: Calculating the initial percentage of H3K4me3 sites.")
print("-" * 50)

# Given values
p_t_fraction = 0.1104  # Proportion of H3K4me3 at time t=10 hours
rate = 0.10             # Decay rate per hour
time = 10               # Time in hours

print("The decay process is modeled by the equation: P(t) = P(0) * e^(-k*t)")
print(f"Given values: P(t) = {p_t_fraction}, k = {rate}, t = {time}")

# Equation setup
print("\nSubstituting the values into the equation:")
print(f"{p_t_fraction} = P(0) * e^(-{rate} * {time})")

# Solve for P(0)
p_0 = p_t_fraction / math.exp(-rate * time)
p_0_percent = p_0 * 100

print("\nSolving for the initial proportion P(0):")
print(f"P(0) = {p_t_fraction} / e^(-{rate * time})")
print(f"P(0) = {p_0:.4f}")
print(f"Therefore, the initial percentage of trimethylated sites is {p_0_percent:.2f}%.\n")


# Part 2: Determine the impact on gene expression
print("Part 2: Calculating the impact on target gene expression.")
print("-" * 50)

# Given values
initial_expression = 200  # RPKM
proportion_decrease = 0.10  # 10% decrease

print("A linear relationship is assumed between H3K4me3 proportion and gene expression.")
print(f"The initial expression is {initial_expression} RPKM, corresponding to the initial proportion of {p_0_percent:.2f}%.")
print(f"The proportion of H3K4me3 sites decreases by {proportion_decrease * 100}%.")

# Calculate the new expression level
# Due to the linear relationship, a 10% decrease in proportion causes a 10% decrease in expression.
new_expression = initial_expression * (1 - proportion_decrease)

print("\nCalculating the new expression level:")
print("New Expression = Initial Expression * (1 - Relative Decrease in Proportion)")
print(f"New Expression = {initial_expression} * (1 - {proportion_decrease})")
print(f"New Expression = {initial_expression} * {1 - proportion_decrease}")
print(f"Therefore, the new gene expression level is {new_expression:.2f} RPKM.")

<<<The initial percentage of trimethylated sites is 30.00%. The new gene expression level is 180.00 RPKM.>>>