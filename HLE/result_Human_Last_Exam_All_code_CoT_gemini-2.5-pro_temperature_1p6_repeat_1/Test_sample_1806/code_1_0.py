def solve_stabilizer_problem():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers are S1=Z1*Z2, S2=Z2*Z3, S3=Z3*Z4.
    """

    # Define stabilizers and logical states using string representations
    # 'Z' for Pauli-Z, 'I' for Identity
    stabilizers = {
        "S1 = Z1*Z2": "ZZII",
        "S2 = Z2*Z3": "IZZI",
        "S3 = Z3*Z4": "IIZZ"
    }
    logical_states = {
        "|0_L> = |0000>": "0000",
        "|1_L> = |1111>": "1111"
    }

    all_checks_passed = True

    # --- Step 1: Check Stabilizer Commutativity ---
    print("--- Step 1: Check Stabilizer Commutativity ---")
    print("The proposed stabilizers are S1=Z1*Z2, S2=Z2*Z3, and S3=Z3*Z4.")
    print("Since the Pauli Z matrix commutes with itself (Z*Z = I) and the Identity matrix (Z*I = Z = I*Z),")
    print("any two stabilizers composed only of Z and I operators will commute with each other.")
    print("For example: [S1, S2] = S1*S2 - S2*S1 = (Z1*Z2)(Z2*Z3) - (Z2*Z3)(Z1*Z2) = Z1*(Z2*Z2)*Z3 - Z2*Z3*Z1*Z2 = Z1*I*Z3 - Z1*Z3 = 0.")
    print("Conclusion: The operators form a valid commuting set.\n")

    # --- Step 2: Check if Logical States are Stabilized ---
    print("--- Step 2: Check if Logical States are Stabilized ---")
    print("We must verify that for every stabilizer S and logical state |psi_L>, S|psi_L> = +1 * |psi_L>.\n")

    # Loop through each logical state
    for state_name, state_str in logical_states.items():
        print(f"--- Checking state {state_name} ---")
        is_fully_stabilized = True

        # Loop through each stabilizer
        for stab_name, stab_str in stabilizers.items():
            total_eigenvalue = 1
            eigenvalue_calc_parts = []
            
            # Calculate the overall eigenvalue by multiplying individual eigenvalues
            for i, pauli_op in enumerate(stab_str):
                if pauli_op == 'Z':
                    qubit_state_char = state_str[i]
                    eigenvalue = 1 if qubit_state_char == '0' else -1
                    eigenvalue_calc_parts.append(str(eigenvalue))
                    total_eigenvalue *= eigenvalue
            
            # Print the detailed calculation for the current stabilizer and state
            print(f"  Action of {stab_name} on {state_name}:")
            # Build the string showing the numbers in the final equation
            final_val_str = " * ".join([f"({v})" for v in eigenvalue_calc_parts])
            print(f"    Calculation of eigenvalue: {final_val_str} = {total_eigenvalue}")

            if total_eigenvalue == 1:
                print(f"    Result: {stab_name} {state_name} = ({total_eigenvalue}) * {state_name}. The state IS stabilized.")
            else:
                print(f"    Result: {stab_name} {state_name} = ({total_eigenvalue}) * {state_name}. The state is NOT stabilized.")
                is_fully_stabilized = False
                all_checks_passed = False
        
        if is_fully_stabilized:
            print(f"-> Conclusion for {state_name}: This state is correctly stabilized by the entire set.\n")
        else:
            print(f"-> Conclusion for {state_name}: This state is NOT correctly stabilized.\n")

    # --- Final Conclusion ---
    print("--- Final Conclusion ---")
    if all_checks_passed:
        print("Yes. Both logical basis states are stabilized by all commuting stabilizers.")
        print("Therefore, the code can be considered a stabilizer code with the given stabilizers.")
    else:
        print("No. Not all logical basis states are stabilized by all the proposed stabilizers.")

solve_stabilizer_problem()
<<<Yes. Both logical basis states are stabilized by all commuting stabilizers.
Therefore, the code can be considered a stabilizer code with the given stabilizers.>>>