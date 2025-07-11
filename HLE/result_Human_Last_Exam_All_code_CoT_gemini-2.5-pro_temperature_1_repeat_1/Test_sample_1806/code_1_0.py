import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    """
    # Define single-qubit states and operators
    q0 = np.array([1, 0], dtype=complex)
    q1 = np.array([0, 1], dtype=complex)
    I = np.identity(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define the logical basis states
    # |0_L> = |0000>
    logical_0 = np.kron(q0, np.kron(q0, np.kron(q0, q0)))
    # |1_L> = |1111>
    logical_1 = np.kron(q1, np.kron(q1, np.kron(q1, q1)))

    # Define the stabilizer generators
    # S1 = Z1 Z2 = Z tensor Z tensor I tensor I
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    # S2 = Z2 Z3 = I tensor Z tensor Z tensor I
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    # S3 = Z3 Z4 = I tensor I tensor Z tensor Z
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

    stabilizers = {'Z1*Z2': S1, 'Z2*Z3': S2, 'Z3*Z4': S3}
    logical_states = {'|0_L>': logical_0, '|1_L>': logical_1}

    print("Step 1: Verify that the logical states are stabilized (have eigenvalue +1).\n")

    all_states_stabilized = True
    for s_name, s_matrix in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")
        for l_name, l_vector in logical_states.items():
            # Apply stabilizer to the logical state
            result_vector = s_matrix @ l_vector

            # Check if the result is the same as the original vector
            is_stabilized = np.allclose(result_vector, l_vector)
            
            # Since the action of Z on |0> or |1> results in +/- 1,
            # the eigenvalue for these states will be +1 or -1.
            eigenvalue = 1 if is_stabilized else -1
            
            # Print the equation showing the result
            print(f"Action: {s_name} on {l_name}")
            print(f"Result: {s_name} {l_name} = {eigenvalue:+d} * {l_name}")

            if not is_stabilized:
                all_states_stabilized = False
                print(f"Verification FAILED: {l_name} is not a +1 eigenvector of {s_name}.\n")
            else:
                print(f"Verification PASSED: {l_name} is stabilized by {s_name}.\n")

    print("\nStep 2: Verify that the stabilizer generators commute.\n")

    # Check commutativity: [S1, S2] = 0?
    commute_12 = np.allclose(S1 @ S2, S2 @ S1)
    print(f"Does Z1*Z2 commute with Z2*Z3? {commute_12}")

    # Check commutativity: [S2, S3] = 0?
    commute_23 = np.allclose(S2 @ S3, S3 @ S2)
    print(f"Does Z2*Z3 commute with Z3*Z4? {commute_23}")
    
    # Check commutativity: [S1, S3] = 0?
    commute_13 = np.allclose(S1 @ S3, S3 @ S1)
    print(f"Does Z1*Z2 commute with Z3*Z4? {commute_13}")

    all_commute = commute_12 and commute_23 and commute_13
    print("\n--- Final Conclusion ---")
    if all_states_stabilized and all_commute:
        print("Yes, this can be considered a stabilizer code with the given stabilizers.")
        print("Both conditions are met:")
        print("1. The logical basis states are +1 eigenvectors of all stabilizers.")
        print("2. The stabilizer generators all commute with each other.")
    else:
        print("No, this cannot be considered a stabilizer code with the given stabilizers.")
        if not all_states_stabilized:
            print("Reason: Not all logical states are stabilized by all generators.")
        if not all_commute:
            print("Reason: The stabilizer generators do not form a commuting set.")

# Run the check
if __name__ == "__main__":
    check_stabilizer_code()