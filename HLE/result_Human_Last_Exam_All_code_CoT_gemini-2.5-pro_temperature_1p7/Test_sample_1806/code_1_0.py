import numpy as np

def main():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    The code is defined by:
    |0_L> = |0000>
    |1_L> = |1111>

    The proposed stabilizers are:
    S1 = Z_1 Z_2
    S2 = Z_2 Z_3
    S3 = Z_3 Z_4
    """

    # --- 1. Define Basic States and Operators ---
    # Computational basis states
    q0 = np.array([[1], [0]])
    q1 = np.array([[0], [1]])

    # Pauli operators
    I = np.eye(2)
    Z = np.array([[1, 0], [0, -1]])

    # --- 2. Define Logical Qubits and Stabilizers ---
    # Logical |0_L> = |0000>
    L0 = np.kron(np.kron(np.kron(q0, q0), q0), q0)
    # Logical |1_L> = |1111>
    L1 = np.kron(np.kron(np.kron(q1, q1), q1), q1)

    # Stabilizer generators
    # S1 = Z_1 Z_2 = Z x Z x I x I
    S1 = np.kron(np.kron(np.kron(Z, Z), I), I)
    # S2 = Z_2 Z_3 = I x Z x Z x I
    S2 = np.kron(np.kron(np.kron(I, Z), Z), I)
    # S3 = Z_3 Z_4 = I x I x Z x Z
    S3 = np.kron(np.kron(np.kron(I, I), Z), Z)

    # Lists for iteration
    stabilizers = {'S1': S1, 'S2': S2, 'S3': S3}
    logical_states = {'|0_L>': L0, '|1_L>': L1}
    
    all_conditions_met = True

    # --- 3. Check Stabilization Condition ---
    print("Step 1: Checking if logical states are stabilized.")
    for s_name, s_op in stabilizers.items():
        for l_name, l_state in logical_states.items():
            # Apply stabilizer to the logical state
            result_state = s_op @ l_state
            
            # Check if S|psi> = |psi>
            is_stabilized = np.allclose(result_state, l_state)
            
            # Print the equation and the result
            print(f"Checking {s_name} {l_name} = +1 * {l_name}: {is_stabilized}")
            
            if not is_stabilized:
                all_conditions_met = False

    # --- 4. Check Commutation Relations ---
    print("\nStep 2: Checking if stabilizer generators commute.")
    stab_items = list(stabilizers.items())
    for i in range(len(stab_items)):
        for j in range(i + 1, len(stab_items)):
            s_i_name, s_i_op = stab_items[i]
            s_j_name, s_j_op = stab_items[j]
            
            # Calculate the commutator: [A, B] = AB - BA
            commutator = (s_i_op @ s_j_op) - (s_j_op @ s_i_op)
            
            # Check if the commutator is the zero matrix
            is_commuting = np.allclose(commutator, np.zeros_like(commutator))
            
            print(f"Checking [{s_i_name}, {s_j_name}] = 0: {is_commuting}")

            if not is_commuting:
                all_conditions_met = False

    # --- 5. Final Conclusion ---
    print("\n--- Conclusion ---")
    if all_conditions_met:
        print("All logical states are stabilized by all generators, and all generators commute.")
        print("Therefore, the code CAN be considered a stabilizer code with the given stabilizers.")
        print("\n<<<Yes>>>")
    else:
        print("Not all conditions were met.")
        print("Therefore, the code CANNOT be considered a stabilizer code with the given stabilizers.")
        print("\n<<<No>>>")

if __name__ == "__main__":
    main()