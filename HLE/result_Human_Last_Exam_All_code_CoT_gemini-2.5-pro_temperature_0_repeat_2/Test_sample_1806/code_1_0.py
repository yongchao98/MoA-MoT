import numpy as np

def main():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    """
    # Define single-qubit states and operators
    ket0 = np.array([1, 0])
    ket1 = np.array([0, 1])
    I = np.identity(2)
    Z = np.array([[1, 0], [0, -1]])

    # Define the 4-qubit logical states
    # |0_L> = |0000>
    ket0L = np.kron(ket0, np.kron(ket0, np.kron(ket0, ket0)))
    # |1_L> = |1111>
    ket1L = np.kron(ket1, np.kron(ket1, np.kron(ket1, ket1)))

    # Define the stabilizer operators
    # S1 = Z_1 Z_2 = Z tensor Z tensor I tensor I
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    # S2 = Z_2 Z_3 = I tensor Z tensor Z tensor I
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    # S3 = Z_3 Z_4 = I tensor I tensor Z tensor Z
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

    stabilizers = {
        "S1 = Z1*Z2": S1,
        "S2 = Z2*Z3": S2,
        "S3 = Z3*Z4": S3
    }

    logical_states = {
        "|0_L>": ket0L,
        "|1_L>": ket1L
    }

    all_stabilized = True

    print("Checking if logical states are stabilized...\n")

    # Check if each stabilizer fixes each logical state
    for s_name, s_op in stabilizers.items():
        for l_name, l_state in logical_states.items():
            # Apply stabilizer to the logical state
            result_state = s_op @ l_state
            
            # The eigenvalue is the first non-zero element of the result
            # divided by the first non-zero element of the original state.
            # For these specific states, we can just look at the first or last element.
            eigenvalue = 0
            if np.any(l_state): # Avoid division by zero for null vectors
                # Find the first non-zero element to calculate the eigenvalue
                idx = np.nonzero(l_state)[0][0]
                eigenvalue = result_state[idx] / l_state[idx]

            print(f"Applying {s_name} to {l_name}:")
            # The final equation is S|psi> = eigenvalue * |psi>
            # We output each part of this equation.
            print(f"{s_name} {l_name} = {eigenvalue.real:+.1f} * {l_name}")

            # Check if the state is an eigenvector with eigenvalue +1
            if not np.allclose(result_state, l_state):
                all_stabilized = False
                print(" -> State is NOT stabilized.\n")
            else:
                print(" -> State is stabilized.\n")

    print("--- Conclusion ---")
    if all_stabilized:
        print("Yes, both logical states |0_L> and |1_L> are stabilized by S1, S2, and S3.")
        print("The code can be considered a stabilizer code with these stabilizers.")
    else:
        print("No, at least one logical state is not stabilized by all operators.")
        print("The code cannot be considered a stabilizer code with these stabilizers.")

if __name__ == "__main__":
    main()