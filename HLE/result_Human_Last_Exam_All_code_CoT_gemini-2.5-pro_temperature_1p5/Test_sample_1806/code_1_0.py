import numpy as np
from functools import reduce

def check_stabilizer_code():
    """
    Checks if the 4-qubit code with |0_L> = |0000> and |1_L> = |1111>
    is a stabilizer code with stabilizers Z1Z2, Z2Z3, and Z3Z4.
    """
    # Define single-qubit states and operators
    q0 = np.array([1, 0], dtype=complex)  # |0>
    q1 = np.array([0, 1], dtype=complex)  # |1>
    op_I = np.identity(2, dtype=complex)   # Identity operator
    op_Z = np.array([[1, 0], [0, -1]], dtype=complex)  # Pauli-Z operator

    # --- 1. Define the logical basis states ---
    # |0_L> = |0000> = |0> x |0> x |0> x |0>
    state_0L = reduce(np.kron, [q0, q0, q0, q0])
    # |1_L> = |1111> = |1> x |1> x |1> x |1>
    state_1L = reduce(np.kron, [q1, q1, q1, q1])

    logical_states = {
        "|0_L>": state_0L,
        "|1_L>": state_1L
    }

    # --- 2. Define the stabilizer generators ---
    # S1 = Z_1 Z_2 = Z x Z x I x I
    S1 = reduce(np.kron, [op_Z, op_Z, op_I, op_I])
    # S2 = Z_2 Z_3 = I x Z x Z x I
    S2 = reduce(np.kron, [op_I, op_Z, op_Z, op_I])
    # S3 = Z_3 Z_4 = I x I x Z x Z
    S3 = reduce(np.kron, [op_I, op_I, op_Z, op_Z])

    stabilizers = {
        "Z1*Z2": S1,
        "Z2*Z3": S2,
        "Z3*Z4": S3
    }
    
    print("Checking if the logical states are stabilized by the generators...\n")
    all_stabilized = True

    # --- 3. Perform the checks ---
    for s_name, S_op in stabilizers.items():
        for psi_name, psi_vec in logical_states.items():
            # Apply stabilizer to the logical state
            result_vec = S_op @ psi_vec
            
            # Check the eigenvalue
            if np.allclose(result_vec, psi_vec):
                eigenvalue = "+1"
            elif np.allclose(result_vec, -psi_vec):
                eigenvalue = "-1"
                all_stabilized = False
            else:
                eigenvalue = "Not an eigenstate"
                all_stabilized = False
            
            # Print the final equation for this check
            print(f"Action: {s_name} {psi_name} = ({eigenvalue}) * {psi_name}")

    # --- 4. Final Conclusion ---
    print("\n--- Conclusion ---")
    if all_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given generators.")
        print("Both logical basis states |0_L> and |1_L> are +1 eigenvectors for all stabilizer generators.")
    else:
        print("No, the code cannot be considered a stabilizer code with the given generators.")
        print("Not all logical basis states are +1 eigenvectors for all stabilizer generators.")

check_stabilizer_code()
<<<Yes>>>