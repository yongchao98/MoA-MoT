import math

# The user can change this value to calculate the result for a different number of qubits.
n = 4

# The problem asks for the number of different destabilizer sets {D_1, ..., D_n}
# for the n-qubit stabilizer set {Z_1, ..., Z_n}.

# The total number of sets is the product of two quantities:
# 1. The number of possible Pauli operator structures for the set {D_1, ..., D_n}.
# 2. The number of ways to assign phase factors to the n generators.

# 1. Number of phase choices:
# Each destabilizer generator D_i can be multiplied by a phase factor from {1, -1, i, -i}.
# Since there are n generators, there are 4^n possible combinations of phases.
num_phases = 4**n

# 2. Number of Pauli structure choices:
# The commutation relations dictate that the choice of the n Pauli strings (ignoring phases)
# corresponds to the choice of an n x n symmetric binary matrix.
# The number of elements defining such a matrix is the sum of elements on and above the diagonal,
# which is n*(n+1)/2.
exponent = n * (n + 1) / 2
# Since each of these elements can be 0 or 1, the total number of such matrices is 2^exponent.
num_pauli_structures = 2**int(exponent)

# The total number of different sets of destabilizers is the product of these two numbers.
total_sets = num_phases * num_pauli_structures

# The final code prints out the components of the calculation as requested.
print(f"For n = {n} qubits:")
print("The total number of different sets of destabilizers is calculated as:")
print("  (Number of phase choices) * (Number of Pauli structure choices)")
print("")
print("The equation is: 4^n * 2^(n*(n+1)/2)")
print(f"-> For n={n}, this becomes: 4^{n} * 2^({n}*({n}+1)/2)")
print(f"-> Which evaluates to: {num_phases} * {num_pauli_structures} = {total_sets}")