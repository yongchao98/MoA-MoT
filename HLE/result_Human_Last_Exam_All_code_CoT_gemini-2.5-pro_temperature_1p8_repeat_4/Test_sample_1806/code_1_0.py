def check_stabilizer_code():
    """
    Checks if a 4-qubit code defined by |0_L> = |0000> and |1_L> = |1111>
    is a stabilizer code with stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.
    """
    # Define the logical states and stabilizer generators
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }
    stabilizers = {
        "S1": ("Z1*Z2", [0, 1]),  # (Name, indices of qubits with Z operator)
        "S2": ("Z2*Z3", [1, 2]),
        "S3": ("Z3*Z4", [2, 3])
    }

    print("--- Verifying Stabilizer Code Conditions ---")
    print("\nCondition 1: Stabilizer generators must commute.")
    print("The stabilizers S1=Z1*Z2, S2=Z2*Z3, S3=Z3*Z4 are products of Pauli Z and Identity matrices.")
    print("Any two such operators always commute. Condition 1 is met.")

    print("\nCondition 2: Logical states must be +1 eigenvectors of all stabilizers.")
    all_stabilized = True

    # Iterate over both logical states
    for state_name, state_bits in logical_states.items():
        print(f"\n--- Checking State {state_name} = |{state_bits}> ---")
        
        # Iterate over the three stabilizers
        for stab_symbol, (stab_ops, stab_indices) in stabilizers.items():
            final_eigenvalue = 1
            eigenvalue_factors = []
            
            # The action of Z on |0> is +1, on |1> is -1.
            # We calculate the total eigenvalue by multiplying the eigenvalues from each
            # qubit the stabilizer acts upon.
            for i in stab_indices:
                bit = state_bits[i]
                if bit == '0':
                    eigenvalue = 1
                    eigenvalue_factors.append(f"({eigenvalue})")
                else:  # bit == '1'
                    eigenvalue = -1
                    eigenvalue_factors.append(f"({eigenvalue})")
                final_eigenvalue *= eigenvalue

            # Create the descriptive output string showing the calculation
            factors_str = " * ".join(eigenvalue_factors)
            
            print(f"Applying {stab_symbol} = {stab_ops}:")
            
            # The final equation requires printing each number
            # This line shows the breakdown of the calculation with each number
            print(f"  Calculation: {stab_ops} |{state_bits}> gives eigenvalues from the Z operators.")
            print(f"  Resulting eigenvalues multiplied: {factors_str} = {final_eigenvalue}")

            # Summarize with the final equation
            equation = f"{stab_symbol} {state_name} = {final_eigenvalue} * {state_name}"
            print(f"  Final Equation: {equation}")
            
            if final_eigenvalue != 1:
                all_stabilized = False
                print("  Conclusion: NOT stabilized.")
            else:
                print("  Conclusion: Stabilized.")

    print("\n--- Final Conclusion ---")
    if all_stabilized:
        print("Yes, all logical basis states are +1 eigenvectors of all stabilizer generators.")
        print("The code CAN be considered a stabilizer code with the given stabilizers.")
    else:
        # This case will not be reached for this specific problem
        print("No, not all logical basis states are +1 eigenvectors of all stabilizer generators.")
        print("The code CANNOT be considered a stabilizer code with the given stabilizers.")

check_stabilizer_code()