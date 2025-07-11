def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers are S1=Z1Z2, S2=Z2Z3, S3=Z3Z4.
    """
    stabilizers = {
        "S1=Z1*Z2": "ZZII",
        "S2=Z2*Z3": "IZZI",
        "S3=Z3*Z4": "IIZZ"
    }
    logical_states = {
        "|0_L> = |0000>": "0000",
        "|1_L> = |1111>": "1111"
    }

    all_stabilized = True

    print("--- Verifying the Stabilizer Conditions ---")

    for state_name, state_str in logical_states.items():
        print(f"\n--- Checking state {state_name} ---")
        for stab_name, stab_str in stabilizers.items():
            phase = 1
            phase_factors = []

            # Calculate the overall phase by applying the stabilizer
            for i in range(4):
                op = stab_str[i]
                qubit_state = state_str[i]
                current_phase = 1
                if op == 'Z' and qubit_state == '1':
                    current_phase = -1
                
                phase *= current_phase
                phase_factors.append(str(current_phase))

            # Print the detailed breakdown of the calculation
            print(f"Applying {stab_name} to {state_name}:")
            equation = f"{stab_name} {state_name} = ({' * '.join(phase_factors)}) * {state_name}"
            print(f"  Calculation: {equation}")
            
            final_result = f"             = {phase} * {state_name}"
            print(final_result)
            
            # Check if the state is stabilized
            if phase == 1:
                print("  Result: The state is stabilized.\n")
            else:
                print(f"  Result: The state is NOT stabilized (Eigenvalue is {phase}).\n")
                all_stabilized = False

    print("--- Final Conclusion ---")
    if all_stabilized:
        print("Yes, both logical basis states |0_L> and |1_L> are stabilized by all three operators.")
        print("Therefore, this code can be considered a stabilizer code with stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.")
    else:
        print("No, at least one logical basis state was not stabilized by one of the operators.")
        print("Therefore, this is NOT a stabilizer code with the given stabilizers.")

if __name__ == "__main__":
    check_stabilizer_code()