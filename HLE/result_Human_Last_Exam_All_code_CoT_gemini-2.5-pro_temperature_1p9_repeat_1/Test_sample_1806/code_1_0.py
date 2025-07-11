import numpy as np

def run_stabilizer_check():
    """
    Checks if a 4-qubit code defined by |0_L> = |0000> and |1_L> = |1111>
    is a stabilizer code with generators Z1Z2, Z2Z3, and Z3Z4.
    """
    # Define basic Pauli matrices and identity
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    I = np.array([[1, 0], [0, 1]], dtype=complex)

    # Define logical basis states |0_L> = |0000> and |1_L> = |1111>
    # |0000> is a vector with 1 at index 0 (binary 0000)
    zero_L = np.zeros(16, dtype=complex)
    zero_L[0] = 1
    # |1111> is a vector with 1 at index 15 (binary 1111)
    one_L = np.zeros(16, dtype=complex)
    one_L[15] = 1

    # Define the stabilizer generators using tensor products (np.kron)
    # S1 = Z1*Z2 = Z tensor Z tensor I tensor I
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))
    # S2 = Z2*Z3 = I tensor Z tensor Z tensor I
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))
    # S3 = Z3*Z4 = I tensor I tensor Z tensor Z
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))

    stabilizers = {'Z1*Z2': S1, 'Z2*Z3': S2, 'Z3*Z4': S3}
    logical_states = {'|0_L>': zero_L, '|1_L>': one_L}

    print("To determine if the code is a stabilizer code with the given stabilizers,")
    print("we must check if the logical basis states are +1 eigenvectors of all stabilizer generators.")
    print("That is, for each stabilizer S and logical state |psi_L>, we must verify S|psi_L> = +1 * |psi_L>.\n")

    all_stabilized = True

    # Iterate through stabilizers and logical states and perform the check
    for s_name, s_matrix in stabilizers.items():
        for l_name, l_vector in logical_states.items():
            print(f"Checking: {s_name} * {l_name}")
            
            # Calculate S |psi_L>
            result_vector = s_matrix @ l_vector
            
            # Check if the result is equal to the original state.
            # np.allclose handles potential floating point inaccuracies.
            is_stabilized = np.allclose(result_vector, l_vector)
            
            # The eigenvalue is +1 if the state is stabilized, -1 otherwise
            # (since stabilizers square to I, eigenvalues are +/-1)
            eigenvalue = 1 if is_stabilized else -1

            # Set the overall flag to False if any check fails
            if not is_stabilized:
                all_stabilized = False
            
            print(f"  Result of multiplication: {s_name}{l_name} = {eigenvalue} * {l_name}")

        print("-" * 30)

    # Print the final conclusion
    print("\nFinal Conclusion:")
    if all_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("All logical basis states are +1 eigenstates of all the stabilizer generators.")
    else:
        print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
        print("At least one logical basis state is not a +1 eigenstate of a stabilizer generator.")

if __name__ == '__main__':
    run_stabilizer_check()