import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code with given stabilizers.
    Logical states: |0_L> = |0000>, |1_L> = |1111>
    Stabilizers: S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4
    """
    # --- Step 1: Define single-qubit states and operators ---
    # Computational basis states
    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)

    # Pauli operators
    I = np.identity(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # --- Step 2: Define logical states and stabilizers for the 4-qubit system ---
    # Logical states using Kronecker product
    ket0L = np.kron(np.kron(np.kron(ket0, ket0), ket0), ket0)
    ket1L = np.kron(np.kron(np.kron(ket1, ket1), ket1), ket1)

    # Stabilizer operators using Kronecker product
    S1 = np.kron(np.kron(np.kron(Z, Z), I), I)  # Z1 * Z2
    S2 = np.kron(np.kron(np.kron(I, Z), Z), I)  # Z2 * Z3
    S3 = np.kron(np.kron(np.kron(I, I), Z), Z)  # Z3 * Z4
    
    stabilizers = {'S1': S1, 'S2': S2, 'S3': S3}
    stabilizer_names = list(stabilizers.keys())

    print("Checking if the code defined by |0_L> = |0000> and |1_L> = |1111>")
    print("can be a stabilizer code with S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4.\n")

    # --- Step 3: Check if stabilizers commute ---
    print("--- Condition 1: Stabilizers must commute ---")
    commute_flag = True
    for i in range(len(stabilizer_names)):
        for j in range(i + 1, len(stabilizer_names)):
            name1 = stabilizer_names[i]
            name2 = stabilizer_names[j]
            op1 = stabilizers[name1]
            op2 = stabilizers[name2]
            
            # Commutator [A, B] = A*B - B*A
            commutator = op1 @ op2 - op2 @ op1
            
            # Check if the commutator is the zero matrix
            is_zero = np.allclose(commutator, np.zeros_like(commutator))
            print(f"[{name1}, {name2}] = 0 : {is_zero}")
            if not is_zero:
                commute_flag = False

    if not commute_flag:
        print("\nConclusion: Not all stabilizers commute. This cannot be a valid stabilizer group.")
        return

    # --- Step 4: Check if logical states are stabilized ---
    print("\n--- Condition 2: Logical states must be stabilized (Eigenvalue = +1) ---")
    
    # Check for |0_L>
    print("\nApplying stabilizers to |0_L> = |0000>:")
    stabilized_0L = True
    for name, S in stabilizers.items():
        result_vec = S @ ket0L
        # Applying Z to |0> gives eigenvalue +1. So Z*Z on |00> gives (+1)*(+1) = 1.
        # Equation: S |0L> = c |0L>. We check if c is 1.
        is_stabilized = np.allclose(result_vec, ket0L)
        print(f"Equation: {name} |0_L> = (+1) * |0_L> : {is_stabilized}")
        if not is_stabilized:
            stabilized_0L = False
    
    # Check for |1_L>
    print("\nApplying stabilizers to |1_L> = |1111>:")
    stabilized_1L = True
    for name, S in stabilizers.items():
        result_vec = S @ ket1L
        # Applying Z to |1> gives eigenvalue -1. So Z*Z on |11> gives (-1)*(-1) = 1.
        # Equation: S |1L> = c |1L>. We check if c is 1.
        is_stabilized = np.allclose(result_vec, ket1L)
        eigenvalue = -1 if np.allclose(result_vec, -ket1L) else 1
        print(f"Equation: {name} |1_L> = ({eigenvalue}) * ({eigenvalue}) * |1_L> = (+1) * |1_L> : {is_stabilized}")
        if not is_stabilized:
            stabilized_1L = False
            
    # --- Step 5: Final Conclusion ---
    print("\n--- Final Conclusion ---")
    if commute_flag and stabilized_0L and stabilized_1L:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("Both conditions are met:")
        print("1. All stabilizers commute with each other.")
        print("2. Both logical basis states, |0000> and |1111>, are stabilized by all stabilizers (have eigenvalue +1).")
        result = "Yes"
    else:
        print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
        result = "No"

    # Suppress the numpy array in the final answer per instruction, just return the text
    return result

# Execute the check
check_stabilizer_code()
