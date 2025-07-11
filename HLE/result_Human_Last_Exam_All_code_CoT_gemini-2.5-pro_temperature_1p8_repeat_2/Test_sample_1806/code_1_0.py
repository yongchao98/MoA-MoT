import numpy as np

def run_stabilizer_check():
    """
    Checks if a 4-qubit code is stabilized by a given set of stabilizers.
    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers are S1=Z1Z2, S2=Z2Z3, S3=Z3Z4.
    """
    # Define Pauli Z and Identity I matrices
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define computational basis states |0> and |1>
    q0 = np.array([1, 0], dtype=complex)
    q1 = np.array([0, 1], dtype=complex)

    # Define the 4-qubit logical basis states using tensor products (np.kron)
    # |0_L> = |0> x |0> x |0> x |0>
    logical_0 = np.kron(q0, np.kron(q0, np.kron(q0, q0)))
    # |1_L> = |1> x |1> x |1> x |1>
    logical_1 = np.kron(q1, np.kron(q1, np.kron(q1, q1)))

    # Define the 4-qubit stabilizer operators using tensor products
    # S1 = Z_1 Z_2 = Z x Z x I x I
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    # S2 = Z_2 Z_3 = I x Z x Z x I
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    # S3 = Z_3 Z_4 = I x I x Z x Z
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

    stabilizers = {'Z1*Z2': S1, 'Z2*Z3': S2, 'Z3*Z4': S3}
    logical_states = {'|0_L>': logical_0, '|1_L>': logical_1}

    print("Checking if logical states are +1 eigenstates of the stabilizer generators.")
    print("-" * 65)

    is_a_stabilizer_code = True

    # Iterate over each stabilizer and logical state to perform the check
    for s_name, S_op in stabilizers.items():
        print(f"Testing stabilizer {s_name}:")
        for l_name, l_state in logical_states.items():
            # Apply the stabilizer operator to the logical state vector
            result_state = np.dot(S_op, l_state)
            
            # The states are computational basis states, so the eigenvalue is a simple sign.
            # S |L> = eigenvalue * |L>. eigenvalue = 1 in all cases here.
            eigenvalue = 1
            
            # Use np.allclose to compare floating-point arrays for equality
            is_stabilized = np.allclose(result_state, l_state)

            print(f"  - Action: {s_name} {l_name} = {eigenvalue:+d} * {l_name}")

            if not is_stabilized:
                print(f"    Result: State {l_name} is NOT stabilized by {s_name}.")
                is_a_stabilizer_code = False
            else:
                print(f"    Result: State {l_name} IS stabilized by {s_name}.")
        print("-" * 65)

    print("\nFinal Conclusion:")
    if is_a_stabilizer_code:
        print("Yes, this is a valid stabilizer code.")
        print("Both logical basis states are simultaneous +1 eigenstates of all commuting stabilizer generators.")
    else:
        print("No, this is not a valid stabilizer code with the given generators.")

if __name__ == '__main__':
    run_stabilizer_check()