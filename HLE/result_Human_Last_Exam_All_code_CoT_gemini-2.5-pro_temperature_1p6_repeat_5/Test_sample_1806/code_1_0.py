import numpy as np

def aprint(arr, precision=3, suppress_small=True):
    """Helper function to pretty-print numpy arrays."""
    with np.printoptions(precision=precision, suppress=suppress_small,
                        formatter={'complex_kind': lambda x: f"{x.real:.{precision}f}{x.imag:+.{precision}f}j"}):
        print(arr)

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with the given stabilizers.
    Logical states: |0_L> = |0000>, |1_L> = |1111>
    Stabilizers: S1=Z1Z2, S2=Z2Z3, S3=Z3Z4
    """
    # --- Step 1: Define basic quantum objects ---
    # Single-qubit states
    ket0 = np.array([1, 0])
    ket1 = np.array([0, 1])

    # Single-qubit operators
    I = np.identity(2)
    Z = np.array([[1, 0], [0, -1]])

    # --- Step 2: Define the logical states and stabilizers for the 4-qubit code ---
    # Logical states
    ket0L = np.kron(np.kron(np.kron(ket0, ket0), ket0), ket0)
    ket1L = np.kron(np.kron(np.kron(ket1, ket1), ket1), ket1)

    # Stabilizer operators
    S1 = np.kron(np.kron(np.kron(Z, Z), I), I)
    S2 = np.kron(np.kron(np.kron(I, Z), Z), I)
    S3 = np.kron(np.kron(np.kron(I, I), Z), Z)

    stabilizers = {'S1': S1, 'S2': S2, 'S3': S3}
    stab_names = list(stabilizers.keys())
    stab_pauli_names = {'S1': 'Z1*Z2', 'S2': 'Z2*Z3', 'S3': 'Z3*Z4'}

    all_conditions_met = True

    # --- Step 3: Check for commutativity ---
    print("--- Verifying Stabilizer Commutativity ---")
    is_commutative = True
    for i in range(len(stab_names)):
        for j in range(i + 1, len(stab_names)):
            name1, name2 = stab_names[i], stab_names[j]
            op1, op2 = stabilizers[name1], stabilizers[name2]
            
            # Commutator C = A*B - B*A
            commutator = op1 @ op2 - op2 @ op1
            
            # Check if the commutator is the zero matrix
            if not np.allclose(commutator, np.zeros_like(commutator)):
                is_commutative = False
                all_conditions_met = False
                print(f"[{name1}, {name2}] != 0. They do NOT commute.")
            else:
                print(f"[{name1}, {name2}] = 0. They commute.")

    if not is_commutative:
        print("\nStabilizer set is not Abelian. This cannot be a stabilizer code.")
    else:
        print("\nAll stabilizers commute. Condition 1 is met.")

    # --- Step 4: Check if logical states are stabilized ---
    print("\n--- Verifying Code Space Stabilization ---")
    is_stabilized = True

    logical_states = {'|0_L>': ket0L, '|1_L>': ket1L}
    logical_state_strings = {'|0_L>': '|0000>', '|1_L>': '|1111>'}

    for stab_name, stabilizer_op in stabilizers.items():
        for state_name, state_vec in logical_states.items():
            # Apply stabilizer to the logical state
            result_vec = stabilizer_op @ state_vec
            
            # Check if result is state_vec (i.e., eigenvalue is +1)
            if not np.allclose(result_vec, state_vec):
                is_stabilized = False
                all_conditions_met = False
                print(f"FAILED: {stab_pauli_names[stab_name]} {logical_state_strings[state_name]} != {logical_state_strings[state_name]}")
            else:
                # Eigenvalue is guaranteed to be real for Hermitian operators
                eigenvalue = np.vdot(state_vec, result_vec).real
                print(f"Check: {stab_pauli_names[stab_name]} {logical_state_strings[state_name]} = ({eigenvalue:+.1f}) {logical_state_strings[state_name]}")

    if not is_stabilized:
        print("\nNot all logical states are stabilized with eigenvalue +1. Condition 2 is NOT met.")
    else:
        print("\nAll logical states are +1 eigenstates of all stabilizers. Condition 2 is met.")
    
    # --- Step 5: Final Conclusion ---
    print("\n--- Conclusion ---")
    if all_conditions_met:
        print("Both conditions are met. The code CAN be considered a stabilizer code with the given stabilizers.")
        final_answer = "Yes"
    else:
        print("One or more conditions failed. The code CANNOT be considered a stabilizer code with the given stabilizers.")
        final_answer = "No"

    return final_answer


# Run the check and print the final answer in the specified format
final_answer = check_stabilizer_code()
print(f"\n<<< {final_answer} >>>")
