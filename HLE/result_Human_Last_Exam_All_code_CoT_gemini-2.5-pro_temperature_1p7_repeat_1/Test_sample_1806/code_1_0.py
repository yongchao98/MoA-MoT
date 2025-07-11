import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    """
    # Define single-qubit states and Pauli matrices
    q0 = np.array([[1], [0]], dtype=complex)  # |0>
    q1 = np.array([[0], [1]], dtype=complex)  # |1>
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # --- Step 1: Define the logical states and stabilizers ---
    
    # Logical states |0_L> = |0000> and |1_L> = |1111>
    L0 = np.kron(np.kron(q0, q0), np.kron(q0, q0))
    L1 = np.kron(np.kron(q1, q1), np.kron(q1, q1))
    
    # Stabilizers S1 = Z_1*Z_2, S2 = Z_2*Z_3, S3 = Z_3*Z_4
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))
    
    stabilizers = {"S1": S1, "S2": S2, "S3": S3}
    all_conditions_met = True

    # --- Step 2: Check if stabilizers commute ---
    print("--- Verifying Stabilizer Commutation ---")
    s_names = list(stabilizers.keys())
    for i in range(len(s_names)):
        for j in range(i + 1, len(s_names)):
            name1, name2 = s_names[i], s_names[j]
            op1, op2 = stabilizers[name1], stabilizers[name2]
            
            # Commutator [S_i, S_j] = S_i*S_j - S_j*S_i
            commutator = op1 @ op2 - op2 @ op1
            
            # Check if the commutator is the zero matrix
            is_commuting = np.allclose(commutator, np.zeros_like(commutator))
            print(f"Does [{name1}, {name2}] = 0?  Result: {is_commuting}")
            if not is_commuting:
                all_conditions_met = False

    # --- Step 3: Check if logical states are stabilized ---
    logical_states = {"|0_L>": L0, "|1_L>": L1}
    for state_name, state_vec in logical_states.items():
        print(f"\n--- Verifying Stabilization of {state_name} ---")
        for s_name, s_op in stabilizers.items():
            # Apply stabilizer to the state: S|psi>
            stabilized_state = s_op @ state_vec
            
            # Check if S|psi> = |psi>
            is_stabilized = np.allclose(stabilized_state, state_vec)
            print(f"Does {s_name}{state_name} = {state_name}? Result: {is_stabilized}")
            if not is_stabilized:
                all_conditions_met = False
    
    # --- Step 4: Final Conclusion ---
    print("\n--- Conclusion ---")
    if all_conditions_met:
        print("Yes, all conditions are met. The code can be considered a stabilizer code with the given stabilizers.")
    else:
        print("No, not all conditions are met. This is not a valid stabilizer code for the given states and stabilizers.")
    
    return "Yes" if all_conditions_met else "No"

# Run the check and print the final answer in the required format
final_answer = check_stabilizer_code()
print(f"\n<<< {final_answer} >>>")