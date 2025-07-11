def check_stabilizer_code():
    """
    Checks if a 4-qubit code can be described by the given stabilizers by applying
    them to the logical basis states and verifying the eigenvalue is +1.
    """
    
    # Define logical states as strings for easy inspection
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }

    # Define stabilizers by the indices (0-based) where Z operators act
    stabilizers = {
        "Z_1 Z_2": [0, 1],
        "Z_2 Z_3": [1, 2],
        "Z_3 Z_4": [2, 3]
    }

    is_valid_code = True

    print("--- Verifying Stabilizer Conditions ---")
    print("A state |psi> is stabilized by S if S|psi> = +1 * |psi>.\n")

    # --- Check |0_L> ---
    state_name_0 = "|0_L>"
    state_str_0 = logical_states[state_name_0]
    print(f"1. Checking the logical state {state_name_0} = |{state_str_0}>:")
    print("   Action of Z on |0> is: Z|0> = +1 * |0>\n")
    
    for stab_name, z_indices in stabilizers.items():
        # For a state like |0000>, any Z operator results in a +1 phase.
        # The product of any number of +1s is still +1.
        eigenvalue = 1
        print(f"   {stab_name} |{state_str_0}> = {eigenvalue:+} |{state_str_0}>")
    print("-" * 35)

    # --- Check |1_L> ---
    state_name_1 = "|1_L>"
    state_str_1 = logical_states[state_name_1]
    print(f"\n2. Checking the logical state {state_name_1} = |{state_str_1}>:")
    print("   Action of Z on |1> is: Z|1> = -1 * |1>\n")

    for stab_name, z_indices in stabilizers.items():
        # For a state like |1111>, each Z operator introduces a -1 phase.
        # Since each stabilizer has two Z operators, the total phase is (-1)*(-1) = +1.
        phase = 1
        factors = []
        for _ in z_indices:
            phase *= -1
            factors.append("(-1)")
        
        eigenvalue = phase
        factors_str = " * ".join(factors)
        
        # Print the full equation showing the calculation of the eigenvalue
        print(f"   {stab_name} |{state_str_1}> = {factors_str} |{state_str_1}> = {eigenvalue:+} |{state_str_1}>")

        if eigenvalue != 1:
            is_valid_code = False
    print("-" * 35)

    # --- Final Conclusion ---
    print("\n--- Conclusion ---")
    if is_valid_code:
        print("Both logical basis states, |0_L> and |1_L>, are stabilized by all three generators.")
        print("Therefore, the code can be considered a stabilizer code with the given stabilizers.")
        final_answer = "Yes"
    else:
        # This part will not be reached for this specific problem
        print("At least one logical basis state was not stabilized by a generator.")
        print("Therefore, this is not a valid stabilizer code with the given stabilizers.")
        final_answer = "No"

    return final_answer

# Execute the check and print the final answer in the required format
final_result = check_stabilizer_code()
print(f"<<<{final_result}>>>")
