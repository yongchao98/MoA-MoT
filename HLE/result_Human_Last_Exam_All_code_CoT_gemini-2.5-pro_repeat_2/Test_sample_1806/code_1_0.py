import numpy as np

def main():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers to check are S1=Z1*Z2, S2=Z2*Z3, S3=Z3*Z4.
    """
    # Define single-qubit states and operators
    I = np.identity(2)
    Z = np.array([[1, 0], [0, -1]])
    q0 = np.array([1, 0])
    q1 = np.array([0, 1])

    # Define the 4-qubit logical states
    L0 = np.kron(q0, np.kron(q0, np.kron(q0, q0)))
    L1 = np.kron(q1, np.kron(q1, np.kron(q1, q1)))
    
    logical_states = {
        "|0_L>": L0,
        "|1_L>": L1
    }

    # Define the 4-qubit stabilizer generators
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))
    
    stabilizers = {
        "Z1*Z2": S1,
        "Z2*Z3": S2,
        "Z3*Z4": S3
    }

    print("--- Verifying if logical states are stabilized ---")
    all_stabilized = True
    for s_name, s_op in stabilizers.items():
        for l_name, l_state in logical_states.items():
            # Apply the stabilizer to the logical state
            result_state = s_op @ l_state
            
            # The eigenvalue lambda is given by <L|S|L>
            # np.vdot is the conjugate dot product
            eigenvalue = np.vdot(l_state, result_state).real
            
            print(f"Applying {s_name} to {l_name}:")
            print(f"Result: {s_name} {l_name} = {eigenvalue:+.1f} * {l_name}")
            
            # Check if the eigenvalue is +1
            if not np.isclose(eigenvalue, 1.0):
                all_stabilized = False
                print(" -> This state is NOT stabilized (eigenvalue is not +1).")
            else:
                 print(" -> This state is stabilized.")
            print("-" * 20)

    print("\n--- Verifying if stabilizer generators commute ---")
    all_commute = True
    s_items = list(stabilizers.items())
    for i in range(len(s_items)):
        for j in range(i + 1, len(s_items)):
            s1_name, s1_op = s_items[i]
            s2_name, s2_op = s_items[j]
            
            # Commutator [A, B] = A*B - B*A
            commutator = (s1_op @ s2_op) - (s2_op @ s1_op)
            
            # Check if the commutator is the zero matrix
            is_zero = np.allclose(commutator, np.zeros_like(commutator))
            print(f"Checking commutator [{s1_name}, {s2_name}]:")
            if is_zero:
                print(f" -> Result: [{s1_name}, {s2_name}] = 0. They commute.")
            else:
                print(f" -> Result: [{s1_name}, {s2_name}] != 0. They DO NOT commute.")
                all_commute = False
            print("-" * 20)
            
    print("\n--- Final Conclusion ---")
    if all_stabilized and all_commute:
        print("Yes, the code can be considered a stabilizer code with the given generators.")
        print("All logical basis states are +1 eigenstates of all generators, and all generators commute.")
    else:
        print("No, the code cannot be considered a stabilizer code with the given generators.")
        if not all_stabilized:
            print("Reason: Not all logical basis states are stabilized by all generators.")
        if not all_commute:
            print("Reason: The stabilizer generators do not all commute with each other.")

if __name__ == "__main__":
    main()