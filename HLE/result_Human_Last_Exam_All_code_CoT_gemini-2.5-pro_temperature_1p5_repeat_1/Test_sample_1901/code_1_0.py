# The design strategy is to start with a promising candidate and tune it.
# Based on step F, we start with n=5 for the alkyl chain C(n)H(2n+1).
n = 5

# Calculate the number of hydrogen atoms in the alkyl chain
num_hydrogens = 2 * n + 1

# Print the general structure and the specific candidate details
print("General molecular structure: C(n)H(2n+1)-Ph-CN")
print(f"Selected candidate based on n = {n}:")

print("\n--- Numbers for the Final Equation ---")
# The "equation" here refers to the parts of the chemical formula.
# We output each number component of the alkyl chain.
print(f"Number of Carbon atoms (n): {n}")
print(f"Number of Hydrogen atoms (2n+1): {num_hydrogens}")

print("\n--- Final Candidate Formula ---")
print(f"The specific molecular formula is: C{n}H{num_hydrogens}-Ph-CN")