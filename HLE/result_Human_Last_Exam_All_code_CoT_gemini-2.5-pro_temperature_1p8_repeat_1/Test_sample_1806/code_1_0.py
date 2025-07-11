import numpy as np

def check_stabilizer(stabilizer_name, stabilizer_op, state_name, state_vec):
    """
    Applies a stabilizer operator to a state vector and checks if it's an
    eigenvector with an eigenvalue of +1. It prints the full result of the check.

    Args:
        stabilizer_name (str): The name of the stabilizer operator (e.g., "S1=Z1Z2").
        stabilizer_op (np.ndarray): The matrix representation of the stabilizer.
        state_name (str): The name of the logical state (e.g., "|0_L>").
        state_vec (np.ndarray): The vector representation of the state.

    Returns:
        bool: True if the state is stabilized (eigenvalue is +1), False otherwise.
    """
    # Apply the stabilizer to the state
    result_vec = stabilizer_op @ state_vec

    # Since the logical states are computational basis states, we can find the eigenvalue
    # by checking the factor by which the state vector was multiplied.
    # We find the index of the non-zero element in the original state vector to get the corresponding element in the result vector.
    try:
        non_zero_idx = np.where(np.abs(state_vec) > 1e-9)[0][0]
        # The eigenvalue is the ratio of the elements at this index.
        eigenvalue = result_vec[non_zero_idx] / state_vec[non_zero_idx]
    except IndexError:
        print(f"Error: State vector {state_name} is a zero vector.")
        return False

    is_stabilized = np.isclose(eigenvalue, 1.0)
    
    # Print the equation showing the result of the operation
    print(f"Applying {stabilizer_name} to {state_name}:")
    print(f"  {stabilizer_name} {state_name} = {eigenvalue.real:+.1f} * {state_name}")
    
    if is_stabilized:
        print(f"  The eigenvalue is +1, so the state is stabilized.")
    else:
        print(f"  The eigenvalue is not +1, so the state is NOT stabilized.")
    print("-" * 50)
    
    return is_stabilized

# --- Main script ---

# Define single-qubit states and operators as numpy arrays
ket0 = np.array([1, 0])
ket1 = np.array([0, 1])
I = np.identity(2)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Construct the 4-qubit logical basis states using the tensor product (np.kron)
# |0_L> = |0000> = |0> ⊗ |0> ⊗ |0> ⊗ |0>
ket0_L = np.kron(np.kron(ket0, ket0), np.kron(ket0, ket0))
# |1_L> = |1111> = |1> ⊗ |1> ⊗ |1> ⊗ |1>
ket1_L = np.kron(np.kron(ket1, ket1), np.kron(ket1, ket1))

# Construct the stabilizer generators using the tensor product
S1 = np.kron(np.kron(Z, Z), np.kron(I, I))  # Z₁Z₂
S2 = np.kron(np.kron(I, Z), np.kron(Z, I))  # Z₂Z₃
S3 = np.kron(np.kron(I, I), np.kron(Z, Z))  # Z₃Z₄

stabilizers = {
    "S1=Z1Z2": S1,
    "S2=Z2Z3": S2,
    "S3=Z3Z4": S3
}
logical_states = {
    "|0_L>": ket0_L,
    "|1_L>": ket1_L
}

# This flag will track if all checks pass
all_stabilized = True

# Perform checks for all combinations of stabilizers and logical states
for s_name, s_op in stabilizers.items():
    for l_name, l_vec in logical_states.items():
        if not check_stabilizer(s_name, s_op, l_name, l_vec):
            all_stabilized = False

# Print the final conclusion based on the results of the checks
print("\n--- Final Conclusion ---")
if all_stabilized:
    print("All logical basis states, |0_L> and |1_L>, are stabilized by all given operators.")
    print("Therefore, the code CAN be considered a stabilizer code with stabilizers Z1Z2, Z2Z3, and Z3Z4.")
else:
    print("At least one logical basis state was not stabilized by one of the operators.")
    print("Therefore, the code CANNOT be considered a stabilizer code with these stabilizers.")