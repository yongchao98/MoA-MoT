import math

# Step 1: Explain the chemical basis for the calculation.
print("The reaction forms a symmetric, dinuclear complex with two identical chiral Ruthenium centers.")
print("This allows us to calculate the number of stereoisomers using a standard formula.")
print("")

# Step 2: Define the number of chiral centers.
n = 2
print(f"Number of equivalent chiral centers (n) = {n}")
print("")

# Step 3: Apply the formula and show the calculation.
# For an even number 'n' of equivalent chiral centers, the formula is:
# Total Isomers = 2^(n-1) + 2^((n/2)-1)
print("Using the formula: Total Isomers = 2^(n-1) + 2^((n/2)-1)")
print("Substituting the value of n:")

term1_exp = n - 1
term2_exp = (n / 2) - 1
term1_val = int(math.pow(2, term1_exp))
term2_val = int(math.pow(2, term2_exp))

print(f"Total Isomers = 2^({n}-1) + 2^(({n}/2)-1)")
print(f"Total Isomers = 2^{int(term1_exp)} + 2^{int(term2_exp)}")
print(f"Total Isomers = {term1_val} + {term2_val}")

# Step 4: Calculate and state the final result.
total_isomers = term1_val + term2_val
print(f"Total Isomers = {total_isomers}")
print("")
print("These three isomers are:")
print("1. The (Δ,Δ) isomer")
print("2. The (Λ,Λ) isomer (enantiomer of the first)")
print("3. The (Δ,Λ) isomer (a meso compound)")
