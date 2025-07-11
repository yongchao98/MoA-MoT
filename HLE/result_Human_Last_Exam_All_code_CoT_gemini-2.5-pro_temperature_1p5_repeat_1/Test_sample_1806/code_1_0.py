import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers are S1=Z1*Z2, S2=Z2*Z3, S3=Z3*Z4.
    """
    # Define single-qubit states and Pauli matrices
    q0 = np.array([1, 0], dtype=complex)
    q1 = np.array([0, 1], dtype=complex)
    I = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define the 4-qubit logical basis states
    state_0L = np.kron(np.kron(np.kron(q0, q0), q0), q0)
    state_1L = np.kron(np.kron(np.kron(q1, q1), q1), q1)
    
    # Define the stabilizer operators
    # S1 = Z_1 * Z_2
    S1 = np.kron(np.kron(np.kron(Z, Z), I), I)
    # S2 = Z_2 * Z_3
    S2 = np.kron(np.kron(np.kron(I, Z), Z), I)
    # S3 = Z_3 * Z_4
    S3 = np.kron(np.kron(np.kron(I, I), Z), Z)

    stabilizers = {
        "S1=Z1*Z2": S1,
        "S2=Z2*Z3": S2,
        "S3=Z3*Z4": S3
    }
    
    all_stabilized = True
    
    print("--- Verifying |0_L> = |0000> ---")
    for name, op in stabilizers.items():
        # Apply stabilizer to the state
        result_state = op @ state_0L
        # The eigenvalue is <psi|O|psi> / <psi|psi>
        # Since state is normalized, eigenvalue = <psi|O|psi>
        eigenvalue = np.vdot(state_0L, result_state).round(2).real
        
        print(f"Applying {name} to |0_L>:")
        # Output the equation with the calculated number
        print(f"{name} |0_L> = {eigenvalue} * |0_L>")
        if not np.isclose(eigenvalue, 1.0):
            print("  -> The state is NOT stabilized (eigenvalue is not +1).")
            all_stabilized = False
        else:
            print("  -> The state is stabilized.")
        print("-" * 20)

    print("\n--- Verifying |1_L> = |1111> ---")
    for name, op in stabilizers.items():
        # Apply stabilizer to the state
        result_state = op @ state_1L
        # Calculate eigenvalue
        eigenvalue = np.vdot(state_1L, result_state).round(2).real

        print(f"Applying {name} to |1_L>:")
        # Output the equation with the calculated number
        print(f"{name} |1_L> = {eigenvalue} * |1_L>")
        if not np.isclose(eigenvalue, 1.0):
            print("  -> The state is NOT stabilized (eigenvalue is not +1).")
            all_stabilized = False
        else:
            print("  -> The state is stabilized.")
        print("-" * 20)

    print("\n--- Conclusion ---")
    if all_stabilized:
        print("Both |0_L> and |1_L> are in the +1 eigenspace of all three stabilizers.")
        print("Therefore, the code can be considered a stabilizer code with the given stabilizers.")
    else:
        print("Not all logical basis states are stabilized by all operators.")
        print("Therefore, this is NOT a stabilizer code with the given stabilizers.")

if __name__ == '__main__':
    check_stabilizer_code()
