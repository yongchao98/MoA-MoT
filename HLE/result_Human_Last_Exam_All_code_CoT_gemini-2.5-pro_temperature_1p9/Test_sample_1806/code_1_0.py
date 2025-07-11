import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    Logical states: |0_L> = |0000>, |1_L> = |1111>
    Stabilizers: S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4
    """
    # Define single-qubit states and operators
    q0 = np.array([1, 0])  # |0>
    q1 = np.array([0, 1])  # |1>
    I = np.identity(2)      # Identity matrix
    Z = np.array([[1, 0], [0, -1]])  # Pauli-Z matrix

    # Define the 4-qubit logical basis states using kronecker product
    L0 = np.kron(np.kron(q0, q0), np.kron(q0, q0)) # |0_L> = |0000>
    L1 = np.kron(np.kron(q1, q1), np.kron(q1, q1)) # |1_L> = |1111>

    # Define the stabilizer operators
    # S1 = Z tensor Z tensor I tensor I
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))
    # S2 = I tensor Z tensor Z tensor I
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))
    # S3 = I tensor I tensor Z tensor Z
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))
    
    stabilizers = {"Z1 Z2": S1, "Z2 Z3": S2, "Z3 Z4": S3}
    logical_states = {"|0_L> = |0000>": L0, "|1_L> = |1111>": L1}
    
    print("--- Verifying Condition 1: Commutativity ---")
    all_commute = True
    s_keys = list(stabilizers.keys())
    for i in range(len(s_keys)):
        for j in range(i + 1, len(s_keys)):
            key1, key2 = s_keys[i], s_keys[j]
            op1, op2 = stabilizers[key1], stabilizers[key2]
            commutator = op1 @ op2 - op2 @ op1
            # Check if commutator is a zero matrix
            if np.allclose(commutator, np.zeros_like(commutator)):
                print(f"[{key1}, {key2}] = 0. They commute.")
            else:
                print(f"[{key1}, {key2}] != 0. They DO NOT commute.")
                all_commute = False
    
    if not all_commute:
        print("\nConclusion: Not a valid stabilizer group as operators do not commute.")
        return False, "No"

    print("\n--- Verifying Condition 2: Stabilizer Action on Logical States ---")
    is_stabilized = True
    for s_name, s_op in stabilizers.items():
        for l_name, l_state in logical_states.items():
            print(f"Applying stabilizer {s_name} to state {l_name}:")
            
            # Apply operator
            result_state = s_op @ l_state
            
            # Check if result_state is a multiple of l_state
            # np.allclose handles floating point comparisons
            if np.allclose(result_state, l_state):
                eigenvalue = 1
            elif np.allclose(result_state, -l_state):
                eigenvalue = -1
            else:
                eigenvalue = "Not an eigenstate" # Should not happen here
            
            print(f"  Equation: {s_name} {l_name} = ({eigenvalue:+}) * {l_name}")
            
            if eigenvalue != 1:
                print(f"  Result: The state is NOT a +1 eigenstate.")
                is_stabilized = False
            else:
                print(f"  Result: The state is a +1 eigenstate.")
            print("-" * 20)

    print("\n--- Final Conclusion ---")
    if all_commute and is_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("Both conditions are met:")
        print("1. All stabilizers commute with each other.")
        print("2. Both logical basis states are +1 eigenstates of all stabilizers.")
        return True, "Yes"
    else:
        print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
        return False, "No"

if __name__ == '__main__':
    is_stabilizer, answer = check_stabilizer_code()
    print(f"\n<<< {answer} >>>")