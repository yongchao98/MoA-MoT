import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with the given generators.
    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers are S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4.
    """
    # Define single-qubit states and operators
    q0 = np.array([1, 0], dtype=complex)
    q1 = np.array([0, 1], dtype=complex)
    I = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # --- Construct logical states (4-qubit states) ---
    # |0_L> = |0> x |0> x |0> x |0>
    L0 = np.kron(np.kron(q0, q0), np.kron(q0, q0))
    # |1_L> = |1> x |1> x |1> x |1>
    L1 = np.kron(np.kron(q1, q1), np.kron(q1, q1))

    # --- Construct stabilizer operators (4-qubit operators) ---
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))  # Z1 * Z2
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))  # Z2 * Z3
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))  # Z3 * Z4

    stabilizers = {
        "S1 = Z1*Z2": S1,
        "S2 = Z2*Z3": S2,
        "S3 = Z3*Z4": S3
    }
    
    logical_states = {
        "|0L> = |0000>": L0,
        "|1L> = |1111>": L1
    }
    
    all_states_stabilized = True
    
    # --- Check 1: Do stabilizers fix the logical states? ---
    print("--- Verifying Stabilizer Conditions ---")
    for state_name, state_vec in logical_states.items():
        print(f"\nChecking state {state_name}...")
        is_state_stabilized = True
        for stab_name, stab_op in stabilizers.items():
            # Apply stabilizer to the state
            result_vec = stab_op @ state_vec
            
            # Eigenvalue is +1 if result_vec is the same as state_vec
            if np.allclose(result_vec, state_vec):
                eigenvalue = 1.0
                print(f"Applying stabilizer {stab_name}: {stab_name.split(' ')[0]}{state_name.split(' ')[0]} = {eigenvalue} * {state_name.split(' ')[0]}")
            # Eigenvalue is -1 if result_vec is -state_vec
            elif np.allclose(result_vec, -state_vec):
                eigenvalue = -1.0
                print(f"Applying stabilizer {stab_name}: {stab_name.split(' ')[0]}{state_name.split(' ')[0]} = {eigenvalue} * {state_name.split(' ')[0]}")
                is_state_stabilized = False
            # Otherwise, it's not an eigenvector
            else:
                print(f"Applying stabilizer {stab_name}: {stab_name.split(' ')[0]}{state_name.split(' ')[0]} is NOT an eigenvector.")
                is_state_stabilized = False
        
        if is_state_stabilized:
            print(f"{state_name.split(' ')[0]} is stabilized by all generators.")
        else:
            print(f"{state_name.split(' ')[0]} is NOT stabilized by all generators.")
            all_states_stabilized = False

    # --- Check 2: Do stabilizers commute? ---
    print("\n--- Verifying Commutation Relations ---")
    stab_list = list(stabilizers.items())
    all_commute = True
    for i in range(len(stab_list)):
        for j in range(i + 1, len(stab_list)):
            name1, op1 = stab_list[i]
            name2, op2 = stab_list[j]
            commutator = (op1 @ op2) - (op2 @ op1)
            
            # Check if commutator is a zero matrix
            if not np.allclose(commutator, np.zeros_like(commutator)):
                all_commute = False
            print(f"[{name1.split(' ')[0]}, {name2.split(' ')[0]}] commutes: {np.allclose(commutator, np.zeros_like(commutator))}")

    # --- Final Conclusion ---
    print("\n--- Conclusion ---")
    if all_states_stabilized and all_commute:
        print("Yes, this code can be considered a stabilizer code with the given stabilizers.")
        global_answer = "Yes"
    else:
        print("No, this code cannot be considered a stabilizer code with the given stabilizers.")
        global_answer = "No"
    return global_answer

if __name__ == '__main__':
    answer = check_stabilizer_code()
    print(f"<<<{answer}>>>")