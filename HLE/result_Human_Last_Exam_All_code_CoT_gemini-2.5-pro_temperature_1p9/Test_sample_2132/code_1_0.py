import math

# Step 1: Define the parameters based on the molecular structure
# Product: methyl cyclopentanecarboxylate (C7H12O2)
N = 21  # Total number of atoms (7 C + 12 H + 2 O)

# Step 2: Define the symmetry classes (Ni)
# The molecule is chiral. The only equivalent atoms are the 3 H's of the methyl group.
# So, we have 1 class of 3 atoms, and (21-3)=18 classes of 1 atom.
# The number of atoms in each class are N_i = [3, 1, 1, 1, ..., 1] (18 times 1)
n_i_values = [3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


# Step 3: Calculate the summation term: Σ(Ni * log2(Ni))
# Note that log2(1) = 0, so classes of size 1 don't contribute to the sum.
sum_term = 0
for n_i in n_i_values:
    if n_i > 1:
        sum_term += n_i * math.log2(n_i)

# Step 4: Calculate the full Böttcher Molecular Complexity
# Formula: C = N*log2(N) - Σ(Ni * log2(Ni))
complexity = N * math.log2(N) - sum_term

# Step 5: Print the results, explaining the calculation
print("Calculation of Böttcher Molecular Complexity for methyl cyclopentanecarboxylate")
print("-" * 70)
print("The formula is: C = N * log2(N) - Σ(Ni * log2(Ni))")
print(f"Total number of atoms (N): {N}")
print("Symmetry classes found:")
print(" - 1 class of 3 atoms (the methyl hydrogens)")
print(" - 18 classes of 1 atom (all other unique atoms)")
print("\nBreaking down the calculation:")

# The final equation as requested by the user prompt
n1 = 3
n2 = 1
count_n2 = 18
print(f"The equation with values is: C = {N} * log2({N}) - ({n1}*log2({n1}) + {count_n2} * ({n2}*log2({n2})))")

term1 = N * math.log2(N)
print(f"Term 1: N * log2(N) = {N} * {math.log2(N):.4f} = {term1:.4f}")
print(f"Term 2 (Summation): Σ(Ni * log2(Ni)) = {n1}*log2({n1}) = {sum_term:.4f}")
print(f"The term '{count_n2} * ({n2}*log2({n2}))' is 0, since log2(1)=0.")

print("\nFinal Result:")
print(f"Complexity = {term1:.4f} - {sum_term:.4f}")
print(f"Böttcher Molecular Complexity = {complexity:.4f}")