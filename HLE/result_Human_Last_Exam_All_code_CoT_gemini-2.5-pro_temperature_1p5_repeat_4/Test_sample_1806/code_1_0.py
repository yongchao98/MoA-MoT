import numpy as np

def print_centered(text, width=70):
    """Helper function to print centered text."""
    print(f"\n{' ' * ((width - len(text)) // 2)}{text}{' ' * ((width - len(text)) // 2)}\n")

def run_stabilizer_check():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    """
    # --- Step 1: Define quantum states and operators ---
    
    # Single-qubit states and operators
    q0 = np.array([[1], [0]], dtype=complex)
    q1 = np.array([[0], [1]], dtype=complex)
    I = np.identity(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Logical basis states for the 4-qubit code
    # |0_L> = |0000>
    state_0L = np.kron(q0, np.kron(q0, np.kron(q0, q0)))
    # |1_L> = |1111>
    state_1L = np.kron(q1, np.kron(q1, np.kron(q1, q1)))
    
    logical_states = {
        "|0_L>": state_0L,
        "|1_L>": state_1L
    }

    # Stabilizer generators as 16x16 matrices
    # S1 = Z_1 Z_2 = Z kron Z kron I kron I
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    # S2 = Z_2 Z_3 = I kron Z kron Z kron I
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    # S3 = Z_3 Z_4 = I kron I kron Z kron Z
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

    stabilizers = {
        "S1 = Z1*Z2": S1,
        "S2 = Z2*Z3": S2,
        "S3 = Z3*Z4": S3
    }
    
    # --- Step 2: Check for commutativity among stabilizers ---
    
    print_centered("--- Verifying Stabilizer Properties ---")
    print("A valid set of stabilizer generators must commute with each other.")
    
    commutes = True
    stab_list = list(stabilizers.items())
    
    for i in range(len(stab_list)):
        for j in range(i + 1, len(stab_list)):
            name1, op1 = stab_list[i]
            name2, op2 = stab_list[j]
            
            # Commutator [A, B] = AB - BA
            commutator = np.dot(op1, op2) - np.dot(op2, op1)
            
            # Check if the commutator is the zero matrix
            is_zero = np.allclose(commutator, np.zeros_like(commutator))
            print(f"Checking [{name1.split(' = ')[0]}, {name2.split(' = ')[0]}]: Result is zero matrix -> {is_zero}")
            if not is_zero:
                commutes = False
    
    if commutes:
        print("\nConclusion: All stabilizer generators commute. This is a valid group.")
    else:
        print("\nConclusion: Not all generators commute. This is not a valid stabilizer group.")
        print("\n<<<No>>>")
        return
        
    # --- Step 3: Check if stabilizers fix the logical states ---

    print_centered("--- Verifying Stabilization of Logical States ---")
    print("The logical states |0_L> and |1_L> must be eigenvectors of each stabilizer with an eigenvalue of +1.")
    
    is_stabilized = True
    for s_name, s_op in stabilizers.items():
        for l_name, l_state in logical_states.items():
            # Apply stabilizer to the logical state
            result_state = np.dot(s_op, l_state)
            
            # Find the index of the first non-zero element to calculate eigenvalue
            # For |0000> this is index 0, for |1111> this is index 15
            first_elem_idx = np.argmax(np.abs(l_state))
            
            # Eigenvalue = (S|psi>)_i / (|psi>)_i
            eigenvalue = result_state[first_elem_idx, 0] / l_state[first_elem_idx, 0]
            
            # Check if S|psi> = eigenvalue * |psi> and if eigenvalue is 1
            is_eigenvector = np.allclose(result_state, eigenvalue * l_state)
            
            print(f"Applying {s_name} to {l_name}:")
            # Output the final equation with the calculated eigenvalue
            print(f"  {s_name.split(' = ')[0]} {l_name} = {eigenvalue.real:.0f} * {l_name}")

            if not (is_eigenvector and np.isclose(eigenvalue, 1)):
                is_stabilized = False

    # --- Step 4: Final Conclusion ---
    print_centered("--- Final Conclusion ---")
    if is_stabilized:
        print("Both logical states are stabilized by all generators (eigenvalue = +1).")
        print("The code can be considered a stabilizer code with the given stabilizers.")
        print("\n<<<Yes>>>")
    else:
        print("Not all logical states are stabilized by all generators.")
        print("The code cannot be considered a stabilizer code with the given stabilizers.")
        print("\n<<<No>>>")

if __name__ == '__main__':
    run_stabilizer_check()