import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code defined by |0_L>=|0000> and |1_L>=|1111>
    is a stabilizer code for S1=Z1Z2, S2=Z2Z3, S3=Z3Z4.
    """
    # Define Pauli matrices
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define logical basis states for a 4-qubit system (16-dimensional space)
    # |0_L> = |0000> which corresponds to the first basis vector (index 0)
    ket0L = np.zeros(16, dtype=complex)
    ket0L[0] = 1
    # |1_L> = |1111> which corresponds to the last basis vector (index 15)
    ket1L = np.zeros(16, dtype=complex)
    ket1L[15] = 1

    logical_states = {
        "|0_L>": ket0L,
        "|1_L>": ket1L
    }
    logical_state_reprs = {
        "|0_L>": "|0000>",
        "|1_L>": "|1111>"
    }

    # Define stabilizer operators using Kronecker product
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))

    stabilizers = {
        "S1 = Z1*Z2": S1,
        "S2 = Z2*Z3": S2,
        "S3 = Z3*Z4": S3
    }

    # Flag to track if all conditions are met
    is_stabilizer_code = True

    print("To be a stabilizer code, all logical basis states must be +1 eigenvectors of all stabilizers.")
    print("Checking condition: S |psi_L> = +1 * |psi_L>\n")

    # Iterate through each stabilizer and logical state
    for s_name, s_op in stabilizers.items():
        print(f"--- Checking Stabilizer {s_name} ---")
        for l_name, l_state in logical_states.items():
            l_repr = logical_state_reprs[l_name]
            # Apply the stabilizer to the logical state
            result_state = s_op @ l_state

            # Since the logical states are computational basis states, the resulting state
            # is just the original state multiplied by the corresponding diagonal element
            # of the stabilizer matrix. We can extract this eigenvalue.
            eigenvalue = 0
            if np.allclose(l_state, ket0L):
                eigenvalue = result_state[0] # Eigenvalue for |0000>
            elif np.allclose(l_state, ket1L):
                eigenvalue = result_state[15] # Eigenvalue for |1111>

            # Ensure the eigenvalue is a real number for printing and comparison
            eigenvalue = np.real(eigenvalue)
            
            # Print the equation as requested
            print(f"Action on {l_name} = {l_repr}:")
            print(f"{s_name} {l_name} = {eigenvalue:+.1f} * {l_name}")

            # Check if the eigenvalue is +1
            if not np.isclose(eigenvalue, 1.0):
                print("  -> Condition FAILED: Eigenvalue is not +1.\n")
                is_stabilizer_code = False
            else:
                print("  -> Condition MET: Eigenvalue is +1.\n")

    print("--- Final Conclusion ---")
    if is_stabilizer_code:
        print("Yes, both logical basis states |0_L> and |1_L> are +1 eigenvectors for all proposed operators.")
        print("The number of physical qubits is n=4 and the number of independent stabilizers is g=3.")
        print("This defines a 2^(n-g) = 2^(4-3) = 2-dimensional stabilizer space.")
        print("The codespace spanned by the two orthogonal logical states is also 2-dimensional and lies within this stabilizer space.")
        print("Therefore, the code can be considered a stabilizer code with the given stabilizers.")
    else:
        print("No, at least one logical basis state is not a +1 eigenvector of all proposed operators.")
        print("Therefore, the code cannot be considered a stabilizer code with the given stabilizers.")

if __name__ == '__main__':
    check_stabilizer_code()