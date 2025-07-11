import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code with |0_L> = |0000> and |1_L> = |1111>
    is a stabilizer code with generators Z1Z2, Z2Z3, and Z3Z4.
    """
    # Define single-qubit states and operators
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    q0 = np.array([[1], [0]], dtype=complex) # |0>
    q1 = np.array([[0], [1]], dtype=complex) # |1>

    # Define the 4-qubit logical basis states using Kronecker product
    # |0_L> = |0> x |0> x |0> x |0> = |0000>
    L0 = np.kron(np.kron(np.kron(q0, q0), q0), q0)
    # |1_L> = |1> x |1> x |1> x |1> = |1111>
    L1 = np.kron(np.kron(np.kron(q1, q1), q1), q1)

    # Define the stabilizer generators as 4-qubit operators
    # S1 = Z_1 Z_2 = Z x Z x I x I
    S1 = np.kron(np.kron(np.kron(Z, Z), I), I)
    # S2 = Z_2 Z_3 = I x Z x Z x I
    S2 = np.kron(np.kron(np.kron(I, Z), Z), I)
    # S3 = Z_3 Z_4 = I x I x Z x Z
    S3 = np.kron(np.kron(np.kron(I, I), Z), Z)

    stabilizers = {'Z1 Z2': S1, 'Z2 Z3': S2, 'Z3 Z4': S3}
    logical_states = {'|0_L>': L0, '|1_L>': L1}
    all_stabilized = True

    print("Checking if the code is a stabilizer code with generators S1=Z1Z2, S2=Z2Z3, S3=Z3Z4.")
    print("A state |psi> is stabilized by an operator S if S|psi> = 1 * |psi>.\n")

    # Iterate through each stabilizer and each logical state to check the condition
    for s_name, S_op in stabilizers.items():
        for l_name, L_state in logical_states.items():
            print(f"Applying S = {s_name} to state {l_name}:")

            # Apply the stabilizer to the logical state
            result_state = S_op @ L_state

            # Find the eigenvalue. For these specific states, the eigenvalue is the first
            # non-zero element of the result state divided by the first non-zero
            # element of the original state.
            eigenvalue = 1.0
            if np.allclose(result_state, L_state):
                eigenvalue = 1.0
            elif np.allclose(result_state, -L_state):
                eigenvalue = -1.0
            else:
                # This case means the state is not an eigenvector
                eigenvalue = "N/A"

            # Print the equation form: S |psi> = lambda |psi>
            print(f"Equation: {s_name} {l_name} = {eigenvalue} * {l_name}")

            # Check for stabilization (eigenvalue must be +1)
            if np.allclose(result_state, L_state):
                print(f"Outcome: {l_name} is stabilized by {s_name}.\n")
            else:
                print(f"Outcome: {l_name} is NOT stabilized by {s_name}.\n")
                all_stabilized = False

    # Print the final conclusion
    print("--------------------------------------------------")
    if all_stabilized:
        print("Final Conclusion: Yes, since all logical basis states are stabilized by all the generators,")
        print("the code defined by |0_L> = |0000> and |1_L> = |1111> can be considered a stabilizer code")
        print("with the stabilizers Z1 Z2, Z2 Z3, and Z3 Z4.")
    else:
        print("Final Conclusion: No, since not all logical basis states are stabilized by all the generators,")
        print("this cannot be considered a stabilizer code with the given stabilizers.")

if __name__ == '__main__':
    check_stabilizer_code()