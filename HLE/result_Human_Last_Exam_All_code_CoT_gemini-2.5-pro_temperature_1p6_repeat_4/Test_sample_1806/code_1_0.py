def check_stabilizer():
    """
    Checks if the logical states |0_L> = |0000> and |1_L> = |1111>
    are stabilized by S1=Z1Z2, S2=Z2Z3, and S3=Z3Z4.
    """
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }
    
    stabilizers = {
        "S1 = Z1*Z2": [0, 1],  # Qubits the stabilizer acts on (0-indexed)
        "S2 = Z2*Z3": [1, 2],
        "S3 = Z3*Z4": [2, 3]
    }

    print("Checking if the logical states are +1 eigenstates of the stabilizers.\n")
    print("Recall: Z|0> = (+1)|0> and Z|1> = (-1)|1>\n")

    all_conditions_met = True

    for s_name, s_qubits in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")
        for l_name, l_state_str in logical_states.items():
            print(f"Applying {s_name} to {l_name} = |{l_state_str}>:")

            # Determine eigenvalues from Z acting on |0> or |1>
            eigen_q1 = 1 if l_state_str[s_qubits[0]] == '0' else -1
            eigen_q2 = 1 if l_state_str[s_qubits[1]] == '0' else -1
            
            # The total eigenvalue is the product of individual eigenvalues
            total_eigenvalue = eigen_q1 * eigen_q2
            
            # Display the calculation
            op_str = f"(Z{s_qubits[0]+1}|{l_state_str[s_qubits[0]]}>)(Z{s_qubits[1]+1}|{l_state_str[s_qubits[1]]}>)"
            result_str = f"({eigen_q1})*({eigen_q2}) |{l_state_str}>"
            
            print(f"{s_name} |{l_state_str}> = {op_str} (acting on the rest trivially)")
            print(f"= {result_str} = {total_eigenvalue}|{l_state_str}>")
            
            if total_eigenvalue == 1:
                print(f"Result: {l_name} is stabilized by {s_name}.\n")
            else:
                print(f"Result: {l_name} is NOT stabilized by {s_name}.\n")
                all_conditions_met = False
    
    print("--- Conclusion ---")
    if all_conditions_met:
        print("Yes, all logical states are +1 eigenstates of all commuting stabilizers.")
        print("Therefore, this can be considered a stabilizer code with the given stabilizers.")
    else:
        print("No, not all logical states are stabilized.")
        print("Therefore, this cannot be considered a stabilizer code with the given stabilizers.")

check_stabilizer()