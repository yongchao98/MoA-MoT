import numpy as np

# Define single-qubit states and operators
q0 = np.array([1, 0], dtype=complex)  # |0>
q1 = np.array([0, 1], dtype=complex)  # |1>
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I = np.identity(2, dtype=complex)

# Define logical states using tensor (Kronecker) product
# |0_L> = |0000>
zero_L = np.kron(np.kron(q0, q0), np.kron(q0, q0))
# |1_L> = |1111>
one_L = np.kron(np.kron(q1, q1), np.kron(q1, q1))

# Define stabilizer generators as matrices
S1 = np.kron(np.kron(Z, Z), np.kron(I, I))  # S1 = Z1*Z2
S2 = np.kron(np.kron(I, Z), np.kron(Z, I))  # S2 = Z2*Z3
S3 = np.kron(np.kron(I, I), np.kron(Z, Z))  # S3 = Z3*Z4

stabilizers = {
    "Z1*Z2": S1,
    "Z2*Z3": S2,
    "Z3*Z4": S3
}

all_stabilized = True

# --- Check for |0_L> ---
print("Checking if |0_L> = |0000> is stabilized:")
for name, S in stabilizers.items():
    result_vector = S @ zero_L
    # Find the eigenvalue by dividing the first non-zero element of the result by the original
    eigenvalue = result_vector[0] / zero_L[0]
    print(f"Action of {name} on |0_L>: {name}|0000> = {eigenvalue.real:.1f} * |0000>")
    if not np.allclose(result_vector, zero_L):
        all_stabilized = False
print("-" * 30)

# --- Check for |1_L> ---
print("Checking if |1_L> = |1111> is stabilized:")
for name, S in stabilizers.items():
    result_vector = S @ one_L
    # Find the eigenvalue by dividing the last non-zero element of the result by the original
    eigenvalue = result_vector[-1] / one_L[-1]
    print(f"Action of {name} on |1_L>: {name}|1111> = {eigenvalue.real:.1f} * |1111>")
    if not np.allclose(result_vector, one_L):
        all_stabilized = False
print("-" * 30)


# --- Final Conclusion ---
print("\nConclusion:")
if all_stabilized:
    print("Yes, the code can be considered a stabilizer code with the given generators.")
    print("Both |0_L> and |1_L> are eigenvectors with an eigenvalue of +1 for all stabilizer operations.")
else:
    print("No, the code cannot be considered a stabilizer code with the given generators.")
