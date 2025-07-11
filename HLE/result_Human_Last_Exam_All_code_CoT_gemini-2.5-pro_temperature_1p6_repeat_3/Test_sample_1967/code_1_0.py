# Set the number of qubits for the system
n = 3

print(f"For an n={n} qubit system with stabilizers S_i = Z_i:")
print("-" * 50)

# Step 1: Count choices for the Pauli part of the destabilizers

# The diagonal elements M_ii of the destabilizer matrix must be X or Y to anti-commute with S_i = Z_i.
# There are n such diagonal elements.
# Number of choices = 2 * 2 * ... * 2 (n times) = 2^n
diagonal_choices = 2**n
print(f"1. Choices for diagonal Pauli operators (M_ii in {{X, Y}}):")
print(f"   The formula is 2^n.")
print(f"   For n={n}, this is 2^{n} = {diagonal_choices}\n")

# The off-diagonal elements M_ij (i!=j) must be I or Z.
# Additionally, for [D_i, D_j]=0 to hold, (M_ij, M_ji) must be either (I, I) or (Z, Z).
# This means for each pair of off-diagonal elements (M_ij, M_ji), there are 2 choices.
# There are n*(n-1)/2 such pairs.
# Number of choices = 2^(n*(n-1)/2)
num_off_diagonal_pairs = n * (n - 1) // 2
off_diagonal_choices = 2**num_off_diagonal_pairs
print(f"2. Choices for off-diagonal Pauli operators (M_ij in {{I, Z}} with M_ij=M_ji mapping):")
print(f"   The formula is 2^(n*(n-1)/2).")
print(f"   For n={n}, this is 2^({num_off_diagonal_pairs}) = {off_diagonal_choices}\n")

# The total number of unique Pauli "skeletons" is the product of the two counts above.
pauli_choices = diagonal_choices * off_diagonal_choices
print(f"3. Total number of distinct Pauli structures:")
print(f"   This is the product of the choices above.")
print(f"   {diagonal_choices} * {off_diagonal_choices} = {pauli_choices}\n")

# Step 2: Count choices for the phases

# Each of the n destabilizers D_i can have a global phase from {+1, -1, +i, -i}.
# Number of choices = 4 * 4 * ... * 4 (n times) = 4^n
phase_choices = 4**n
print(f"4. Choices for global phases (4 for each of n destabilizers):")
print(f"   The formula is 4^n.")
print(f"   For n={n}, this is 4^{n} = {phase_choices}\n")

# Step 3: Calculate the total number of different destabilizer sets

# The total is the number of Pauli skeletons times the number of phase combinations.
total_choices = pauli_choices * phase_choices
print("-" * 50)
print("The total number of different sets of destabilizers is the final product:")
print(f"Final Equation: ({diagonal_choices} * {off_diagonal_choices}) * {phase_choices}")
print(f"Total = {pauli_choices} * {phase_choices} = {total_choices}")
