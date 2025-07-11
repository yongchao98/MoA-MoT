def solve():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    """
    # Define logical states as lists of bits (0 or 1)
    logical_states = {
        "|0_L> = |0000>": [0, 0, 0, 0],
        "|1_L> = |1111>": [1, 1, 1, 1]
    }

    # Define stabilizer generators as lists of Pauli operators ('Z' or 'I')
    # S1 = Z1 Z2 I3 I4
    # S2 = I1 Z2 Z3 I4
    # S3 = I1 I2 Z3 Z4
    stabilizers = {
        "S1 = Z1*Z2": ['Z', 'Z', 'I', 'I'],
        "S2 = Z2*Z3": ['I', 'Z', 'Z', 'I'],
        "S3 = Z3*Z4": ['I', 'I', 'Z', 'Z']
    }

    print("To be a stabilizer code, all logical states must be stabilized by all generators.")
    print("This means S|psi> = +1 * |psi> for every generator S and logical state |psi>.\n")

    is_stabilizer_code = True

    # Iterate over each stabilizer
    for s_name, s_ops in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")
        
        # Iterate over each logical state
        for l_name, l_bits in logical_states.items():
            total_eigenvalue = 1
            eigenvalue_factors = []
            
            # Iterate over each qubit (from 1 to 4)
            for i in range(4):
                op = s_ops[i]
                bit = l_bits[i]
                eigenvalue = 1
                
                # Z|0> = +1|0>, Z|1> = -1|1>
                if op == 'Z':
                    if bit == 1:
                        eigenvalue = -1
                
                # I|psi> = +1|psi>, so we don't need to do anything for 'I'
                
                total_eigenvalue *= eigenvalue
                # Add each number to the equation
                eigenvalue_factors.append(str(eigenvalue))
            
            # Format the equation string
            equation_str = " * ".join([f"({x})" for x in eigenvalue_factors])
            
            print(f"Applying {s_name} to {l_name}:")
            print(f"  Eigenvalue contributions from (qubit 1, 2, 3, 4) respectively:")
            # Here we print each number in the final equation
            print(f"  Final eigenvalue = {equation_str} = {total_eigenvalue}")
            
            # Check if the state is stabilized (i.e., eigenvalue is +1)
            if total_eigenvalue == 1:
                print("  Result: STABILIZED")
            else:
                is_stabilizer_code = False
                print(f"  Result: NOT STABILIZED (Eigenvalue is {total_eigenvalue})")
            
        print("")

    # Final conclusion
    print("--- Conclusion ---")
    if is_stabilizer_code:
        print("Yes. As shown above, both logical basis states |0_L> and |1_L> have an eigenvalue of +1 for all three generators.")
        print("Therefore, this code can be considered a stabilizer code with the stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.")
    else:
        print("No. At least one logical basis state is not stabilized by one of the operators (eigenvalue is not +1).")
        print("Therefore, this code CANNOT be considered a stabilizer code with the given generators.")

solve()