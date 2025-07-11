import numpy as np

def verify_stabilizer_code():
    """
    Verifies if a 4-qubit code defined by |0L>=|0000> and |1L>=|1111>
    is a stabilizer code with stabilizers Z1Z2, Z2Z3, Z3Z4.
    """
    # --- Step 1: Define single-qubit states and operators ---
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)

    # --- Step 2: Construct the logical basis states using Kronecker products ---
    # |0L> = |0> x |0> x |0> x |0>
    ket0_L = np.kron(ket0, np.kron(ket0, np.kron(ket0, ket0)))
    # |1L> = |1> x |1> x |1> x |1>
    ket1_L = np.kron(ket1, np.kron(ket1, np.kron(ket1, ket1)))

    logical_states = {'|0L>': ket0_L, '|1L>': ket1_L}
    
    print("Logical states defined as:")
    print("|0L> = |0000>")
    print("|1L> = |1111>\n")

    # --- Step 3: Construct the stabilizer operators using Kronecker products ---
    # S1 = Z_1 * Z_2 * I_3 * I_4
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    # S2 = I_1 * Z_2 * Z_3 * I_4
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    # S3 = I_1 * I_2 * Z_3 * Z_4
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

    stabilizers = {'Z1*Z2': S1, 'Z2*Z3': S2, 'Z3*Z4': S3}
    
    print("Proposed stabilizers:")
    for name in stabilizers:
        print(f"S = {name}")
    print("\nVerifying that the logical states are +1 eigenstates of the stabilizers...")
    
    # --- Step 4: Verify the stabilizer conditions and print results ---
    all_stabilized = True
    for l_name, l_state in logical_states.items():
        print(f"\n--- Checking logical state {l_name} ---")
        for s_name, s_op in stabilizers.items():
            # Apply stabilizer S to logical state |psi_L>
            result_state = s_op @ l_state
            
            # The eigenvalue lambda is given by <psi_L|S|psi_L> since |psi_L> is normalized.
            eigenvalue = np.vdot(l_state, result_state).real
            
            # Print the equation S |psi_L> = lambda |psi_L>
            print(f"({s_name}) {l_name} = {eigenvalue:+.0f} * {l_name}")

            if not np.isclose(eigenvalue, 1.0):
                all_stabilized = False
    
    # --- Step 5: Final Conclusion ---
    print("\n--- Conclusion ---")
    if all_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("Both |0L> and |1L> are simultaneous +1 eigenvectors for all proposed stabilizers,")
        print("and the stabilizers form a commuting group.")
    else:
        print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
        print("At least one logical state is not a +1 eigenvector of all stabilizers.")

if __name__ == '__main__':
    verify_stabilizer_code()