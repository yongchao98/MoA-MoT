def check_stabilizer_code():
    """
    Checks if the 4-qubit repetition code is stabilized by the given operators.
    The script iterates through each stabilizer and each logical state, calculates
    the eigenvalue, and prints the result of the operation.
    """
    stabilizers = {
        "S1 = Z1*Z2": ['Z', 'Z', 'I', 'I'],
        "S2 = Z2*Z3": ['I', 'Z', 'Z', 'I'],
        "S3 = Z3*Z4": ['I', 'I', 'Z', 'Z']
    }
    logical_states = {
        "|0_L> = |0000>": "0000",
        "|1_L> = |1111>": "1111"
    }

    print("To be a stabilizer code, each logical basis state must be an eigenstate with eigenvalue +1 for every stabilizer.")
    print("We check the condition: S |psi_L> = +1 * |psi_L>\n")
    
    all_stabilized = True

    # Iterate through each stabilizer
    for s_name, s_paulis in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")
        
        # Iterate through each logical state
        for l_name, l_state_str in logical_states.items():
            print(f"Applying {s_name} to {l_name}:")
            
            total_eigenvalue = 1
            eigenvalue_nums = []
            
            # Apply the operator to each qubit in the state
            for i in range(4):
                pauli_op = s_paulis[i]
                qubit_state = l_state_str[i]
                
                current_eigenvalue = 1
                if pauli_op == 'Z' and qubit_state == '1':
                    current_eigenvalue = -1
                
                total_eigenvalue *= current_eigenvalue
                eigenvalue_nums.append(current_eigenvalue)

            # Format the numbers for the final equation output
            eigenvalue_strs = [f"({x})" for x in eigenvalue_nums]
            product_str = " * ".join(eigenvalue_strs)
            
            # Print the step-by-step calculation
            print(f"  The total eigenvalue is the product of the individual eigenvalues:")
            print(f"  {product_str} = {total_eigenvalue}")

            if total_eigenvalue == 1:
                print(f"  Final Equation: {s_name} {l_name} = {total_eigenvalue} * {l_name}")
                print(f"  Result: {l_name} is stabilized by {s_name}.\n")
            else:
                print(f"  Final Equation: {s_name} {l_name} = {total_eigenvalue} * {l_name}")
                print(f"  Result: {l_name} is NOT stabilized by {s_name}.\n")
                all_stabilized = False
                
    print("--- Overall Conclusion ---")
    if all_stabilized:
        print("Yes. Both |0_L> and |1_L> are eigenstates with eigenvalue +1 for all three operators.")
        print("Therefore, this is a valid stabilizer code with the given stabilizers.")
    else:
        print("No. At least one logical basis state is not an eigenstate with eigenvalue +1 for all operators.")

# Execute the check
check_stabilizer_code()