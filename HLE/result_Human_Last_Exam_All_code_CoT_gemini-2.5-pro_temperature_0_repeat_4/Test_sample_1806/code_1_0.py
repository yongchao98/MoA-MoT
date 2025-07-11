import numpy as np

def main():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    """
    # Define single-qubit states and operators
    q0 = np.array([1, 0], dtype=complex)
    q1 = np.array([0, 1], dtype=complex)
    I = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define the logical basis states
    # |0_L> = |0000>
    L0 = np.kron(q0, np.kron(q0, np.kron(q0, q0)))
    # |1_L> = |1111>
    L1 = np.kron(q1, np.kron(q1, np.kron(q1, q1)))

    # Define the stabilizer generators
    # S1 = Z1 Z2
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    # S2 = Z2 Z3
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    # S3 = Z3 Z4
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))

    stabilizers = {
        "S1 = Z1*Z2": S1,
        "S2 = Z2*Z3": S2,
        "S3 = Z3*Z4": S3
    }
    
    logical_states = {
        "|0_L> = |0000>": L0,
        "|1_L> = |1111>": L1
    }

    print("Step 1: Checking if stabilizer generators commute.")
    s_keys = list(stabilizers.keys())
    commute_ok = True
    for i in range(len(s_keys)):
        for j in range(i + 1, len(s_keys)):
            name1 = s_keys[i].split(" = ")[0]
            name2 = s_keys[j].split(" = ")[0]
            op1 = stabilizers[s_keys[i]]
            op2 = stabilizers[s_keys[j]]
            
            # Commutator [A, B] = AB - BA
            commutator = op1 @ op2 - op2 @ op1
            if not np.allclose(commutator, np.zeros_like(commutator)):
                print(f"[{name1}, {name2}] != 0. They do not commute.")
                commute_ok = False
            else:
                print(f"[{name1}, {name2}] = 0. They commute.")
    
    if not commute_ok:
        print("\nSince not all stabilizers commute, this is not a valid stabilizer group.")
        print("\nConclusion: No, the code cannot be considered a stabilizer code with these generators.")
        return

    print("\nAll stabilizer generators commute. This is a valid stabilizer group.")
    print("\nStep 2: Checking if logical states are stabilized.")
    
    all_stabilized = True
    
    for s_name, s_op in stabilizers.items():
        for l_name, l_state in logical_states.items():
            print("-" * 30)
            print(f"Applying {s_name.split(' = ')[0]} to {l_name.split(' = ')[0]}:")
            
            # Apply stabilizer to the logical state
            result_state = s_op @ l_state
            
            # Check if the resulting state is the same as the original
            # np.allclose handles potential floating point inaccuracies
            is_stabilized = np.allclose(result_state, l_state)
            
            # Determine the eigenvalue
            eigenvalue = 0
            if is_stabilized:
                eigenvalue = 1.0
            elif np.allclose(result_state, -l_state):
                eigenvalue = -1.0
                all_stabilized = False
            else:
                # This case shouldn't happen for these specific inputs
                eigenvalue = "Not an eigenvector"
                all_stabilized = False

            print(f"Equation: {s_name.split(' = ')[0]} {l_name.split(' = ')[0]} = {eigenvalue} * {l_name.split(' = ')[0]}")

            if is_stabilized:
                print(f"Result: {l_name.split(' = ')[0]} is stabilized by {s_name.split(' = ')[0]}.")
            else:
                print(f"Result: {l_name.split(' = ')[0]} is NOT stabilized by {s_name.split(' = ')[0]}.")
                all_stabilized = False

    print("=" * 40)
    print("\nFinal Conclusion:")
    if all_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("Both logical states are stabilized by all generators, and the generators commute.")
    else:
        print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
        print("At least one logical state was not stabilized by all generators.")

if __name__ == "__main__":
    main()
<<<Yes, the code can be considered a stabilizer code with the given stabilizers.>>>