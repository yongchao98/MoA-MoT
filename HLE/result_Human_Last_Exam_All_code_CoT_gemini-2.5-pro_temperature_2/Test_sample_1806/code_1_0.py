import numpy as np

def check_pauli_commutation(pauli_str1, pauli_str2):
    """
    Checks if two Pauli strings commute. They commute if they anticommute on an even number of qubits.
    For this problem, since we only have Z and I operators, they will always commute.
    """
    anticommuting_qubits = 0
    anticommute_pairs = {('X', 'Y'), ('Y', 'X'), ('X', 'Z'), ('Z', 'X'), ('Y', 'Z'), ('Z', 'Y')}

    for p1, p2 in zip(pauli_str1, pauli_str2):
        if (p1, p2) in anticommute_pairs:
            anticommuting_qubits += 1

    return anticommuting_qubits % 2 == 0

def run_checks():
    """
    Runs the verification checks for the quantum error-correcting code.
    """
    # Define the stabilizers and logical states
    stabilizers = {
        "S1=Z1*Z2": "ZZII",
        "S2=Z2*Z3": "IZZI",
        "S3=Z3*Z4": "IIZZ"
    }
    logical_states = {
        "|0_L>": "0000",
        "|1_L>": "1111"
    }

    # --- Step 1: Check if stabilizers commute ---
    print("Step 1: Checking if stabilizers commute.")
    all_commute = True
    stabilizer_items = list(stabilizers.items())

    for i in range(len(stabilizer_items)):
        for j in range(i + 1, len(stabilizer_items)):
            name1, op1 = stabilizer_items[i]
            name2, op2 = stabilizer_items[j]
            commute = check_pauli_commutation(op1, op2)
            print(f"  Checking [{name1.split('=')[0]}, {name2.split('=')[0]}]: {'Commute' if commute else 'Anticommute'}")
            if not commute:
                all_commute = False

    if all_commute:
        print("Result: All stabilizer generators commute with each other.\n")
    else:
        print("Result: Not all stabilizers commute. This cannot be a valid stabilizer code.\n")

    # --- Step 2: Check if logical states are stabilized ---
    print("Step 2: Checking if logical states are stabilized (eigenvalue must be +1).")
    all_stabilized = True

    for s_name, s_op_str in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")

        for l_name, l_state_str in logical_states.items():
            print(f"  Action on {l_name} = |{l_state_str}>:")
            
            eigenvalues = []
            final_eigenvalue = 1
            op_str_parts = []
            
            # Action of Z on |0> is +1, on |1> is -1
            for k, (pauli_char, bit_char) in enumerate(zip(s_op_str, l_state_str)):
                if pauli_char == 'Z':
                    op_str_parts.append(f"Z{k+1}")
                    if bit_char == '0':
                        eigenvalues.append(1)
                    else: # bit_char == '1'
                        eigenvalue = -1
                        eigenvalues.append(eigenvalue)
                        final_eigenvalue *= eigenvalue

            # Constructing and printing the equation
            op_display = "".join(op_str_parts)
            eigenvalue_display = " * ".join([f"({v})" for v in eigenvalues])
            print(f"    {op_display}|{l_state_str}> = ({eigenvalue_display})|{l_state_str}> = {final_eigenvalue}|{l_state_str}>")
            
            if final_eigenvalue != 1:
                all_stabilized = False
                print(f"    Result: {l_name} is NOT stabilized.\n")
            else:
                print(f"    Result: {l_name} is stabilized.\n")

    # --- Step 3: Final Conclusion ---
    print("--- Final Conclusion ---")
    if all_commute and all_stabilized:
        print("Yes, the code defined by |0_L> = |0000> and |1_L> = |1111> can be considered a stabilizer code with stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.")
    else:
        print("No, this cannot be considered a stabilizer code with the given stabilizers.")

# Execute the verification
run_checks()