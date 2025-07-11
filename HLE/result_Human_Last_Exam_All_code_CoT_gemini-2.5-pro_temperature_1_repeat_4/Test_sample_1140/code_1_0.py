# Plan:
# 1. Define the parameters for the synthesis: peptide length and the efficiency of each amino acid coupling step.
# 2. A 100-amino-acid peptide requires 99 coupling steps.
# 3. The overall theoretical yield is calculated by raising the single-step efficiency to the power of the number of steps.
# 4. Print the equation and the final result to demonstrate the challenge of synthesizing long peptides directly.

# 1. Define parameters
peptide_length = 100
# A high, yet realistic, coupling efficiency for a single step in SPPS
coupling_efficiency = 0.995 

# 2. Calculate the number of coupling steps
num_couplings = peptide_length - 1

# 3. Calculate the overall theoretical yield
# The formula is: Overall Yield = (Coupling Efficiency) ^ (Number of Couplings)
overall_yield = coupling_efficiency ** num_couplings

# 4. Print the explanation, the equation with numbers, and the result.
print(f"Calculating the theoretical yield for a {peptide_length}aa peptide with a {coupling_efficiency*100}% single-step efficiency:")
print("-" * 20)
# This line fulfills the requirement to output each number in the final equation
print(f"Equation: Overall Yield = {coupling_efficiency} ^ {num_couplings}")
print(f"Resulting Overall Yield: {overall_yield:.4f}")
print(f"Resulting Overall Yield as a percentage: {overall_yield * 100:.2f}%")
print("-" * 20)
print("\nThis low theoretical yield demonstrates why multi-fragment approaches like Native Chemical Ligation are preferred for long peptides.")
