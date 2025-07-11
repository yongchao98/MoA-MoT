import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    Code: |0_L> = |0000>, |1_L> = |1111>
    Stabilizers: S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4
    """
    # --- Define single-qubit states and operators ---
    q0 = np.array([[1], [0]], dtype=complex)
    q1 = np.array([[0], [1]], dtype=complex)
    I = np.identity(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # --- Construct logical states using tensor products (np.kron) ---
    L0 = np.kron(np.kron(q0, q0), np.kron(q0, q0))
    L1 = np.kron(np.kron(q1, q1), np.kron(q1, q1))
    
    logical_states = [
        ('|0_L> = |0000>', L0),
        ('|1_L> = |1111>', L1)
    ]

    # --- Construct stabilizer operators using tensor products ---
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))
    
    stabilizers = [
        ('S1 = Z1*Z2', S1),
        ('S2 = Z2*Z3', S2),
        ('S3 = Z3*Z4', S3)
    ]

    print("--- Checking if logical states are stabilized ---")
    all_states_stabilized = True
    for s_name, S in stabilizers:
        for l_name, L in logical_states:
            # Apply stabilizer to the logical state
            result_state = S @ L
            
            # Check if the state is an eigenvector and find the eigenvalue
            eigenvalue = 0
            if np.allclose(result_state, L):
                eigenvalue = 1
            elif np.allclose(result_state, -L):
                eigenvalue = -1
            
            is_stabilized = (eigenvalue == 1)
            if not is_stabilized:
                all_states_stabilized = False

            # Print the equation for this check
            print(f"Checking {s_name} on {l_name}:")
            print(f"Result: {s_name} {l_name} = {eigenvalue} * {l_name}")
            print(f"Is it stabilized (eigenvalue=+1)? {'Yes' if is_stabilized else 'No'}\n")

    print("\n--- Checking if stabilizer generators commute ---")
    all_generators_commute = True
    # Check pairs (S1, S2), (S1, S3), (S2, S3)
    pairs = [
        (stabilizers[0], stabilizers[1]),
        (stabilizers[0], stabilizers[2]),
        (stabilizers[1], stabilizers[2])
    ]
    for (s1_name, S1_op), (s2_name, S2_op) in pairs:
        # Calculate commutator: S1*S2 - S2*S1
        commutator = (S1_op @ S2_op) - (S2_op @ S1_op)
        
        # Check if the commutator is the zero matrix
        commutes = np.allclose(commutator, np.zeros_like(commutator))
        if not commutes:
            all_generators_commute = False
            
        print(f"Checking [{s1_name}, {s2_name}]:")
        print(f"Result: [{s1_name}, {s2_name}] = 0? {'Yes' if commutes else 'No'}\n")

    # --- Final Conclusion ---
    print("--- Final Conclusion ---")
    if all_states_stabilized and all_generators_commute:
        print("Yes, all logical states are stabilized by all generators, and all generators commute.")
        print("The code can be considered a stabilizer code with the given stabilizers.")
        final_answer = "Yes"
    else:
        print("No, one or more conditions were not met.")
        final_answer = "No"
        
    return final_answer

if __name__ == '__main__':
    answer = check_stabilizer_code()
    # The final answer is enclosed in <<< >>> as requested.
    # print(f"\n<<<{answer}>>>") # This would be for the final submission format.

# Execute the function to show the step-by-step analysis.
check_stabilizer_code()