def check_stabilizer_code():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.

    The code is defined by:
    |0_L> = |0000>
    |1_L> = |1111>

    The proposed stabilizers are:
    S1 = Z1 * Z2
    S2 = Z2 * Z3
    S3 = Z3 * Z4
    """

    # Define logical states as strings
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }

    # Define stabilizers as lists of operators ('Z' or 'I' for Identity)
    stabilizers = {
        "S1 = Z1*Z2": ['Z', 'Z', 'I', 'I'],
        "S2 = Z2*Z3": ['I', 'Z', 'Z', 'I'],
        "S3 = Z3*Z4": ['I', 'I', 'Z', 'Z']
    }

    print("--- Verifying the Stabilizer Conditions ---")
    print("A state |psi> is stabilized by an operator S if S|psi> = (+1)|psi>.\n")

    is_stabilizer_code = True

    # Iterate through each logical state
    for state_name, state_bits in logical_states.items():
        print(f"--- Checking State: {state_name} = |{state_bits}> ---\n")
        
        # Iterate through each stabilizer
        for stab_name, stab_ops in stabilizers.items():
            eigenvalue = 1
            
            # Apply the stabilizer to the state
            # The Pauli Z operator flips the sign for a |1> state, eigenvalue -1
            # and leaves the sign for a |0> state, eigenvalue +1
            for i in range(len(state_bits)):
                if stab_ops[i] == 'Z' and state_bits[i] == '1':
                    eigenvalue *= -1
            
            # Print the equation for this check
            print(f"Applying {stab_name}:")
            # The state itself doesn't change because Z is diagonal in the computational basis
            print(f"{stab_name} {state_name} = {stab_name} |{state_bits}> = {eigenvalue:+d} * |{state_bits}>")

            if eigenvalue != 1:
                is_stabilizer_code = False
                print(f"Result: {state_name} is NOT stabilized by {stab_name}.\n")
            else:
                print(f"Result: {state_name} IS stabilized by {stab_name}.\n")

    print("--- Final Conclusion ---")
    if is_stabilizer_code:
        print("Yes, the code can be considered a stabilizer code with the stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.")
        print("This is because both |0_L> and |1_L> are +1 eigenstates of all three operators.")
    else:
        print("No, this cannot be considered a stabilizer code with the given stabilizers.")
        print("At least one logical state was not a +1 eigenstate of one of the operators.")

if __name__ == "__main__":
    check_stabilizer_code()
<<<Yes, the code can be considered a stabilizer code with the stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.>>>