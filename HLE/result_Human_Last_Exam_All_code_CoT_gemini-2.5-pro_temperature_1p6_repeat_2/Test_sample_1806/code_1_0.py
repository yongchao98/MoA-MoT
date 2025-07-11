def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with the given stabilizers.
    The logical states are |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers are S1 = Z1*Z2, S2 = Z2*Z3, and S3 = Z3*Z4.
    """
    # Define the stabilizer generators by the indices of the Z operators
    stabilizers = {
        "Z1Z2": [0, 1],
        "Z2Z3": [1, 2],
        "Z3Z4": [2, 3]
    }

    # Define the logical basis states as strings
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }

    # This flag will be set to False if any state is not stabilized
    is_a_stabilizer_code = True

    print("Checking if the code is stabilized by the given generators...\n")

    # The action of Z on |0> gives eigenvalue +1, and on |1> gives -1.
    def get_z_eigenvalue(qubit_state):
        return 1 if qubit_state == '0' else -1

    # Iterate over both logical states
    for state_name, state_str in logical_states.items():
        print(f"Applying stabilizers to state {state_name} = |{state_str}>:")
        
        # Iterate over the three stabilizer generators
        for stab_name, qubit_indices in stabilizers.items():
            # The total eigenvalue is the product of individual eigenvalues
            total_eigenvalue = 1
            for index in qubit_indices:
                total_eigenvalue *= get_z_eigenvalue(state_str[index])

            # Print the result of the operation in equation form
            # The prompt requires outputting the number in the final equation.
            # Here, total_eigenvalue is the number.
            print(f"{stab_name} |{state_str}> = {total_eigenvalue:+} |{state_str}>")

            # Check if the state is actually stabilized (eigenvalue must be +1)
            if total_eigenvalue != 1:
                is_a_stabilizer_code = False
        print("-" * 30)

    # Print the final conclusion
    if is_a_stabilizer_code:
        print("Conclusion: Yes. Both logical basis states have an eigenvalue of +1 for all stabilizer generators.")
        print("This means the codespace is stabilized, and the code can be considered a stabilizer code with these generators.")
    else:
        print("Conclusion: No. At least one logical basis state is not stabilized by a generator (eigenvalue is not +1).")
        print("Therefore, this cannot be considered a stabilizer code with these generators.")

if __name__ == "__main__":
    check_stabilizer_code()
