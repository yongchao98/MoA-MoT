import numpy as np

# A stabilizer code requires that all logical states |psi_L> are +1 eigenstates
# of the stabilizer generators S_i. We will check this condition: S_i |psi_L> = +1 * |psi_L>.

# --- Define quantum mechanical objects ---

# Pauli Z matrix and 2x2 Identity matrix
Z = np.array([[1, 0], [0, -1]])
I = np.identity(2)

# Computational basis states |0> and |1>
q0 = np.array([1, 0])
q1 = np.array([0, 1])

# Logical basis states |0_L> = |0000> and |1_L> = |1111>
# We use np.kron for the tensor product.
L0 = np.kron(q0, np.kron(q0, np.kron(q0, q0)))
L1 = np.kron(q1, np.kron(q1, np.kron(q1, q1)))

# Stabilizer generators S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4
S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

# --- Perform the verification ---

stabilizers = {
    'S1 = Z1*Z2': S1,
    'S2 = Z2*Z3': S2,
    'S3 = Z3*Z4': S3
}
logical_states = {
    '|0_L>': L0,
    '|1_L>': L1
}

all_stabilized = True

print("Checking if logical states are stabilized...\n")

# Iterate through each stabilizer and each logical state
for s_name, s_matrix in stabilizers.items():
    for l_name, l_vector in logical_states.items():
        # Apply the stabilizer operator to the logical state vector
        result_vector = s_matrix @ l_vector

        # Check if the resulting vector is the same as the original vector.
        # If it is, the eigenvalue is +1.
        if np.allclose(result_vector, l_vector):
            eigenvalue = 1
        else:
            # If not, the state is not stabilized with eigenvalue +1.
            eigenvalue = -1 # In this problem, it would be -1 or some other value
            all_stabilized = False

        print(f"Action of {s_name} on {l_name}:")
        # The prompt requires printing the number in the final equation.
        # The equation is S|psi> = lambda|psi>, so we print lambda.
        print(f"{s_name} {l_name} = {eigenvalue} * {l_name}\n")

# --- Final Conclusion ---
print("--------------------------------------------------")
if all_stabilized:
    print("Conclusion: YES.")
    print("Both logical basis states, |0_L> and |1_L>, are +1 eigenvectors of all three generators.")
    print("Therefore, the code can be considered a stabilizer code with these generators.")
else:
    print("Conclusion: NO.")
    print("At least one logical basis state is not a +1 eigenvector of all generators.")
