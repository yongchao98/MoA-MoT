import math

# Part 1: Calculate the initial percentage of trimethylated sites

print("--- Part 1: Initial H3K4me3 Percentage Calculation ---")

# Given values
p_final = 0.1104  # Final proportion after 10 hours
k = 0.10          # Decay rate of 10% per hour
t = 10            # Time in hours

print("The problem is modeled using the exponential decay formula: P(t) = P_0 * (1 - k)^t")
print(f"We are given P({t}) = {p_final}, k = {k}, and t = {t}.")
print("To find the initial proportion P_0, we rearrange the formula: P_0 = P(t) / (1 - k)^t")
print("\nStep 1: Calculate the decay factor (1 - k)^t")
decay_base = 1 - k
decay_factor = math.pow(decay_base, t)
print(f"(1 - {k})^{t} = ({decay_base})^{t} = {decay_factor:.6f}")

print("\nStep 2: Calculate the initial proportion P_0")
p_initial = p_final / decay_factor
print(f"P_0 = {p_final} / {decay_factor:.6f}")
print(f"P_0 = {p_initial:.6f}")

# Convert to percentage and display the result
p_initial_percent = p_initial * 100
print(f"\nThe initial proportion of trimethylated sites is {p_initial_percent:.2f}%.")

print("\n" + "="*50 + "\n")

# Part 2: Calculate the impact on gene expression

print("--- Part 2: Gene Expression Impact Calculation ---")

# Given values
e_initial = 200.0  # Initial gene expression in RPKM
decrease_fraction = 0.10 # 10% decrease

print("A linear relationship exists between H3K4me3 proportion and gene expression.")
print("Therefore, a 10% decrease in the proportion results in a 10% decrease in expression.")
print("\nStep 1: Define the equation for the new expression level")
print("New Expression = Initial Expression * (1 - Decrease Fraction)")
print(f"Initial Expression = {e_initial} RPKM")
print(f"Decrease Fraction = {decrease_fraction}")

print("\nStep 2: Calculate the new expression level")
e_new = e_initial * (1 - decrease_fraction)
print(f"New Expression = {e_initial} * (1 - {decrease_fraction})")
print(f"New Expression = {e_initial} * {1 - decrease_fraction}")
print(f"New Expression = {e_new:.2f} RPKM")

print(f"\nThe impact is a reduction in target gene expression to {e_new:.2f} RPKM.")

# Final combined answer in the required format
final_answer = f"Initial percentage: {p_initial_percent:.2f}%, New expression: {e_new:.2f} RPKM"
# The prompt asks for one answer, but there are two questions.
# The second question is arguably the main "impact" question. So returning the second numerical answer.
final_numeric_answer = e_new

# Let's present both answers as the prompt implies there might be two results.
# However, the example format <<<9.8>>> suggests a single number. The second question about the "impact" is the final part of the prompt.
# So I will output the final gene expression level.

# print(f"<<<{final_answer}>>>")
# print(f"<<<{final_numeric_answer}>>>")