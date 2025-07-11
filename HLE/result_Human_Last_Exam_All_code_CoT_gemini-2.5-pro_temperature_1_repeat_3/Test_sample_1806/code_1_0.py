def check_stabilizer_code():
    """
    Checks if the code defined by |0_L> = |0000> and |1_L> = |1111>
    is a stabilizer code with generators Z1Z2, Z2Z3, and Z3Z4.
    """
    # Define logical states as strings for easy bit access
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }

    # Define stabilizers by their name and the indices (0-based) of the qubits they act on
    stabilizers = [
        ("Z1 Z2", (0, 1)),
        ("Z2 Z3", (1, 2)),
        ("Z3 Z4", (2, 3))
    ]

    print("To be a valid stabilizer code, the logical states must be +1 eigenvectors of all stabilizers.")
    print("-" * 70)

    all_stabilized = True

    def get_z_action_eigenvalue(bit):
        """Returns the eigenvalue of Z acting on a computational basis state."""
        return 1 if bit == '0' else -1

    # Iterate over each logical state and check if it's stabilized by all generators
    for state_name, state_str in logical_states.items():
        print(f"Checking state {state_name} = |{state_str}>:")
        is_current_state_stabilized = True
        
        for stab_name, indices in stabilizers:
            # Calculate the overall eigenvalue by multiplying eigenvalues from each Z operator
            eigenvalue = 1
            factors = []
            for qubit_index in indices:
                bit = state_str[qubit_index]
                factor = get_z_action_eigenvalue(bit)
                eigenvalue *= factor
                factors.append(str(factor))

            # Format and print the equation showing the calculation
            factor_str = " * ".join([f"({f})" for f in factors])
            print(f"  Action of {stab_name}:")
            print(f"  {stab_name} {state_name} = {factor_str} * {state_name} = {eigenvalue} * {state_name}")
            
            if eigenvalue != 1:
                is_current_state_stabilized = False
                all_stabilized = False
        
        if is_current_state_stabilized:
            print(f"  Result: {state_name} is stabilized by all generators.\n")
        else:
            print(f"  Result: {state_name} is NOT stabilized by all generators.\n")

    print("-" * 70)
    print("Conclusion:")
    if all_stabilized:
        print("Both |0_L> and |1_L> are stabilized by all proposed generators.")
        print("Since the generators also form a valid commuting set, the answer is Yes.")
    else:
        print("At least one logical state was not stabilized by a generator.")
        print("Therefore, the answer is No.")

check_stabilizer_code()