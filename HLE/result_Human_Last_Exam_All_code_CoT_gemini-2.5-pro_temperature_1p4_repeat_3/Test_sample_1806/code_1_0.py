import numpy as np

def verify_stabilizer_code():
    """
    Verifies if a 4-qubit code is a stabilizer code with given stabilizers.
    """
    # Define single-qubit states and matrices
    q0 = np.array([1, 0])
    q1 = np.array([0, 1])
    I = np.identity(2)
    Z = np.array([[1, 0], [0, -1]])

    # Define logical basis states |0_L> = |0000> and |1_L> = |1111>
    L0 = np.kron(np.kron(q0, q0), np.kron(q0, q0))
    L1 = np.kron(np.kron(q1, q1), np.kron(q1, q1))
    
    logical_states = {
        "|0000>": L0,
        "|1111>": L1
    }

    # Define stabilizer generators
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))

    stabilizers = {
        "Z1 Z2": S1,
        "Z2 Z3": S2,
        "Z3 Z4": S3
    }
    
    all_stabilized = True

    # Iterate through each stabilizer and each logical state to check the condition S|psi> = |psi>
    for s_name, S in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")
        for l_name, L_state in logical_states.items():
            # Apply stabilizer to the logical state
            result_state = S @ L_state
            
            # Check if the resulting state is the same as the original state (eigenvalue +1)
            if np.allclose(result_state, L_state):
                eigenvalue = 1
                print(f"Applying {s_name} on {l_name}:")
                # Print the equation as required
                print(f"{s_name} {l_name} = {eigenvalue} * {l_name}")
            else:
                # This part will not be reached if the code is correct
                all_stabilized = False
                print(f"Error: {s_name} does not stabilize {l_name}.")

    print("\n--- Conclusion ---")
    if all_stabilized:
        print("Yes, the code defined by |0_L> = |0000> and |1_L> = |1111> can be considered a stabilizer code with stabilizers Z1 Z2, Z2 Z3, and Z3 Z4.")
    else:
        print("No, this is not a valid stabilizer code for the given logical states.")

if __name__ == "__main__":
    verify_stabilizer_code()