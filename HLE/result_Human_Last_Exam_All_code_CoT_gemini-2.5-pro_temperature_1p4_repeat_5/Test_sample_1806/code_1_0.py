import numpy as np

def run_stabilizer_check():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    The code is defined by |0L> = |0000> and |1L> = |1111>.
    The stabilizers are S1=Z1Z2, S2=Z2Z3, S3=Z3Z4.
    """
    # Define Pauli matrices and Identity
    I = np.identity(2)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define basis kets for a single qubit
    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)

    # Define 4-qubit logical states using Kronecker product
    # |0L> = |0000>
    ket0L = np.kron(ket0, np.kron(ket0, np.kron(ket0, ket0)))
    # |1L> = |1111>
    ket1L = np.kron(ket1, np.kron(ket1, np.kron(ket1, ket1)))

    # Define stabilizer operators as 16x16 matrices
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

    stabilizers = {"S1 = Z1*Z2": S1, "S2 = Z2*Z3": S2, "S3 = Z3*Z4": S3}
    logical_states = {"|0L> = |0000>": ket0L, "|1L> = |1111>": ket1L}

    is_stabilizer_code = True

    # --- Step 1: Check if logical states are stabilized ---
    print("Step 1: Checking if the logical states are stabilized by the operators.\n")
    for s_name, s_op in stabilizers.items():
        s_symbol = s_name.split('=')[0].strip()
        for l_name, l_state in logical_states.items():
            l_symbol = l_name.split('=')[0].strip()
            
            # Apply the stabilizer to the logical state
            result_state = s_op @ l_state
            
            # The eigenvalue is <psi|S|psi> / <psi|psi>. For normalized basis states, it's just <psi|S|psi>.
            # It will be a 1x1 matrix, so we extract the single complex value.
            eigenvalue = (l_state.conj().T @ result_state)[0, 0]

            print(f"Applying {s_symbol} to {l_symbol}:")
            # We output the real part of the eigenvalue, which should be an integer (+1 or -1)
            print(f"Equation: {s_symbol} {l_symbol} = {int(np.real(eigenvalue))} * {l_symbol}")
            
            # Check if the state is unchanged (eigenvalue is +1)
            if not np.allclose(result_state, l_state):
                print(f"Result: {l_symbol} is NOT stabilized by {s_symbol}.\n")
                is_stabilizer_code = False
            else:
                print(f"Result: {l_symbol} is stabilized by {s_symbol}.\n")

    # --- Step 2: Check if stabilizer generators commute ---
    print("Step 2: Checking if the stabilizer generators commute with each other.\n")
    s_list = list(stabilizers.items())
    for i in range(len(s_list)):
        for j in range(i + 1, len(s_list)):
            s_i_name = s_list[i][0].split('=')[0].strip()
            s_j_name = s_list[j][0].split('=')[0].strip()
            s_i_op = s_list[i][1]
            s_j_op = s_list[j][1]
            
            # Calculate the commutator: [A, B] = AB - BA
            commutator = s_i_op @ s_j_op - s_j_op @ s_i_op
            
            print(f"Checking commutator [{s_i_name}, {s_j_name}]:")
            if not np.allclose(commutator, np.zeros_like(commutator)):
                print(f"Result: {s_i_name} and {s_j_name} DO NOT commute.\n")
                is_stabilizer_code = False
            else:
                print(f"Result: {s_i_name} and {s_j_name} commute.\n")

    # --- Step 3: Final Conclusion ---
    print("--- Conclusion ---")
    if is_stabilizer_code:
        print("Yes, this is a valid stabilizer code. Both conditions (stabilization and commutation) are met.")
    else:
        print("No, this is not a valid stabilizer code because one or more conditions failed.")

if __name__ == '__main__':
    run_stabilizer_check()