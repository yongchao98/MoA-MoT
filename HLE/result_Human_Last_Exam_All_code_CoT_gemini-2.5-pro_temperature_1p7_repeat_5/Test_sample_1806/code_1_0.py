import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with the given generators.
    Logical states: |0_L> = |0000>, |1_L> = |1111>
    Stabilizers: S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4
    """
    # Define single-qubit basis states and Pauli matrices (as column vectors/matrices)
    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)
    I = np.identity(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Helper function to create tensor products of operators
    def op_tensor(ops):
        t = np.array([[1.0]]) # Start with a 1x1 identity
        for op in ops:
            t = np.kron(t, op)
        return t

    # --- Step 1: Define the logical states and stabilizers ---
    print("--- Defining States and Stabilizers ---")
    
    # Define logical states
    ket0L = op_tensor([ket0, ket0, ket0, ket0])
    ket1L = op_tensor([ket1, ket1, ket1, ket1])
    print("Logical state |0_L> defined as |0000>")
    print("Logical state |1_L> defined as |1111>\n")

    # Define stabilizer generators
    S1 = op_tensor([Z, Z, I, I])
    S2 = op_tensor([I, Z, Z, I])
    S3 = op_tensor([I, I, Z, Z])
    
    stabilizers = {
        "Z1*Z2": S1,
        "Z2*Z3": S2,
        "Z3*Z4": S3
    }
    print("Stabilizer generators defined: Z1*Z2, Z2*Z3, Z3*Z4\n")
    
    all_conditions_met = True

    # --- Step 2: Check if stabilizer generators commute ---
    print("--- Checking Commutation of Stabilizers ---")
    
    # Check [S1, S2]
    comm_S1_S2 = S1 @ S2 - S2 @ S1
    if not np.allclose(comm_S1_S2, np.zeros((16, 16))):
        all_conditions_met = False
    print(f"[Z1*Z2, Z2*Z3] = 0 : {all_conditions_met}")

    # Check [S1, S3]
    comm_S1_S3 = S1 @ S3 - S3 @ S1
    s1s3_commutes = np.allclose(comm_S1_S3, np.zeros((16, 16)))
    if not s1s3_commutes:
        all_conditions_met = False
    print(f"[Z1*Z2, Z3*Z4] = 0 : {s1s3_commutes}")
        
    # Check [S2, S3]
    comm_S2_S3 = S2 @ S3 - S3 @ S2
    s2s3_commutes = np.allclose(comm_S2_S3, np.zeros((16, 16)))
    if not s2s3_commutes:
        all_conditions_met = False
    print(f"[Z2*Z3, Z3*Z4] = 0 : {s2s3_commutes}\n")
    
    if not all_conditions_met:
        print("Error: Stabilizer generators do not commute.\n")


    # --- Step 3: Check if logical states are stabilized ---
    print("--- Checking if Logical States are Stabilized (S|psi> = c|psi>) ---")
    for name, S in stabilizers.items():
        # Check |0_L> = |0000>
        # Eigenvalue c = <psi|S|psi> since states are normalized
        eigenvalue_0L = (ket0L.conj().T @ S @ ket0L)[0, 0]
        print(f"For {name} acting on |0_L> = |0000>:")
        print(f"The equation is {name}|0000> = {eigenvalue_0L.real:.1f} * |0000>")
        if not np.isclose(eigenvalue_0L, 1.0):
            all_conditions_met = False

        # Check |1_L> = |1111>
        eigenvalue_1L = (ket1L.conj().T @ S @ ket1L)[0, 0]
        print(f"For {name} acting on |1_L> = |1111>:")
        print(f"The equation is {name}|1111> = {eigenvalue_1L.real:.1f} * |1111>")
        if not np.isclose(eigenvalue_1L, 1.0):
            all_conditions_met = False
        print("-" * 30)

    # --- Step 4: Final Conclusion ---
    print("\n--- Final Conclusion ---")
    if all_conditions_met:
        print("Yes, this code can be considered a stabilizer code with the given stabilizers.")
        print("All generators commute and stabilize both logical basis states with an eigenvalue of +1.")
    else:
        print("No, this code cannot be considered a stabilizer code with the given stabilizers.")
        print("One or more of the required conditions were not met.")

if __name__ == '__main__':
    check_stabilizer_code()