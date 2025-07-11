import numpy as np

def get_operator(op_string):
    """Creates a multi-qubit operator matrix from a string description."""
    pauli_matrices = {
        'I': np.array([[1, 0], [0, 1]], dtype=complex),
        'X': np.array([[0, 1], [1, 0]], dtype=complex),
        'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
        'Z': np.array([[1, 0], [0, -1]], dtype=complex),
    }
    
    op_list = [pauli_matrices[p] for p in op_string]
    
    if not op_list:
        return None
        
    # Build the full operator using kronecker products
    full_op = op_list[0]
    for op in op_list[1:]:
        full_op = np.kron(full_op, op)
        
    return full_op

def main():
    """
    Checks if a code can be described by a given set of stabilizers.
    """
    # Define stabilizers
    stabilizer_strings = ["ZZII", "IZZI", "IIZZ"]
    stabilizers = {f"S{i+1}": get_operator(s) for i, s in enumerate(stabilizer_strings)}
    
    # Define logical states |0_L> = |0000> and |1_L> = |1111>
    L0 = np.zeros(16, dtype=complex)
    L0[0] = 1
    L0 = L0.reshape(-1, 1)

    L1 = np.zeros(16, dtype=complex)
    L1[-1] = 1
    L1 = L1.reshape(-1, 1)

    # --- 1. Commutation Check ---
    print("--- 1. Commutation Check ---")
    all_commute = True
    s_names = list(stabilizers.keys())
    for i in range(len(s_names)):
        for j in range(i + 1, len(s_names)):
            name1 = s_names[i]
            name2 = s_names[j]
            S1 = stabilizers[name1]
            S2 = stabilizers[name2]
            
            commutator = S1 @ S2 - S2 @ S1
            
            if not np.allclose(commutator, np.zeros((16, 16))):
                all_commute = False
            
            print(f"[{name1}, {name2}] = 0: {np.allclose(commutator, np.zeros((16, 16)))}")
    print("-" * 30)
    
    # --- 2. Stabilization Check ---
    print("\n--- 2. Stabilization Check ---")
    all_stabilized = True
    for name, S in stabilizers.items():
        # Check |0_L>
        res_L0 = S @ L0
        # The eigenvalue is the dot product of the normalized resulting vector with the original
        eigenvalue_L0 = np.vdot(L0, res_L0) 
        if not np.allclose(res_L0, L0):
            all_stabilized = False
        print(f"Applying {name} on |0_L>: {name} |0_L> = ({np.round(eigenvalue_L0, 2)}) * |0_L>.  State is stabilized: {np.allclose(res_L0, L0)}")
        
        # Check |1_L>
        res_L1 = S @ L1
        eigenvalue_L1 = np.vdot(L1, res_L1)
        if not np.allclose(res_L1, L1):
            all_stabilized = False
        print(f"Applying {name} on |1_L>: {name} |1_L> = ({np.round(eigenvalue_L1, 2)}) * |1_L>.  State is stabilized: {np.allclose(res_L1, L1)}")
        print()
    print("-" * 30)

    # --- 3. Final Conclusion ---
    print("\n--- 3. Final Conclusion ---")
    if all_commute and all_stabilized:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("Reason: All stabilizers commute with each other, and they all stabilize the logical basis states |0_L> and |1_L> with an eigenvalue of +1.")
    else:
        print("No, the code cannot be considered a stabilizer code with the given stabilizers.")
        if not all_commute:
            print("Reason: Not all stabilizers commute.")
        if not all_stabilized:
            print("Reason: Not all logical states are stabilized by all stabilizers.")

if __name__ == "__main__":
    main()
