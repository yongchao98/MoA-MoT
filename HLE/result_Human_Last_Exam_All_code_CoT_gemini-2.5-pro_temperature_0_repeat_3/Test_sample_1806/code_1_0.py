import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code can be described by a given set of stabilizers.
    The code is defined by |0_L> = |0000> and |1_L> = |1111>.
    The stabilizers are S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4.
    """
    # Define single-qubit operators and states
    I = np.identity(2)
    Z = np.array([[1, 0], [0, -1]])
    q0 = np.array([[1], [0]])  # |0>
    q1 = np.array([[0], [1]])  # |1>

    # Define 4-qubit logical states using the Kronecker product
    # |0_L> = |0000>
    L0 = np.kron(q0, np.kron(q0, np.kron(q0, q0)))
    # |1_L> = |1111>
    L1 = np.kron(q1, np.kron(q1, np.kron(q1, q1)))

    # Define 4-qubit stabilizer operators
    S1 = np.kron(Z, np.kron(Z, np.kron(I, I)))  # Z1*Z2
    S2 = np.kron(I, np.kron(Z, np.kron(Z, I)))  # Z2*Z3
    S3 = np.kron(I, np.kron(I, np.kron(Z, Z)))  # Z3*Z4

    stabilizers = {'Z1*Z2': S1, 'Z2*Z3': S2, 'Z3*Z4': S3}
    logical_states = {'|0_L>': L0, '|1_L>': L1}

    all_stabilized = True

    print("--- Checking Stabilization Condition: S|psi> = |psi> ---\n")

    # Check if all logical states are stabilized by all stabilizers
    for s_name, S in stabilizers.items():
        for l_name, L in logical_states.items():
            # Apply stabilizer to the logical state
            result_state = S @ L
            is_stabilized = np.allclose(result_state, L)

            print(f"Applying {s_name} to {l_name}:")
            # Print the symbolic equation showing the eigenvalues
            if l_name == '|0_L>':
                if s_name == 'Z1*Z2':
                    print(f"Equation: {s_name} |0000> = (Z1|0>)(Z2|0>)|00> = (+1)|0> * (+1)|0> * |00> = 1 * |0000>")
                elif s_name == 'Z2*Z3':
                    print(f"Equation: {s_name} |0000> = |0>(Z2|0>)(Z3|0>)|0> = |0> * (+1)|0> * (+1)|0> * |0> = 1 * |0000>")
                elif s_name == 'Z3*Z4':
                    print(f"Equation: {s_name} |0000> = |00>(Z3|0>)(Z4|0>) = |00> * (+1)|0> * (+1)|0> = 1 * |0000>")
            elif l_name == '|1_L>':
                if s_name == 'Z1*Z2':
                    print(f"Equation: {s_name} |1111> = (Z1|1>)(Z2|1>)|11> = (-1)|1> * (-1)|1> * |11> = 1 * |1111>")
                elif s_name == 'Z2*Z3':
                    print(f"Equation: {s_name} |1111> = |1>(Z2|1>)(Z3|1>)|1> = |1> * (-1)|1> * (-1)|1> * |1> = 1 * |1111>")
                elif s_name == 'Z3*Z4':
                    print(f"Equation: {s_name} |1111> = |11>(Z3|1>)(Z4|1>) = |11> * (-1)|1> * (-1)|1> = 1 * |1111>")

            print(f"Result: The state is stabilized. (Verified numerically: {is_stabilized})\n")

            if not is_stabilized:
                all_stabilized = False

    print("\n--- Checking Commutation Condition: [Si, Sj] = 0 ---\n")

    # Check if stabilizers commute
    commute12 = np.allclose(S1 @ S2, S2 @ S1)
    print(f"Checking if Z1*Z2 and Z2*Z3 commute: [S1, S2] = 0? {commute12}")

    commute23 = np.allclose(S2 @ S3, S3 @ S2)
    print(f"Checking if Z2*Z3 and Z3*Z4 commute: [S2, S3] = 0? {commute23}")

    commute13 = np.allclose(S1 @ S3, S3 @ S1)
    print(f"Checking if Z1*Z2 and Z3*Z4 commute: [S1, S3] = 0? {commute13}\n")

    all_commute = commute12 and commute23 and commute13

    print("--- Conclusion ---")
    if all_stabilized and all_commute:
        print("Yes, the code can be considered a stabilizer code with the given stabilizers.")
        print("This is because both logical basis states are stabilized by all generators, and all generators commute with each other.")
        global_answer = "Yes"
    else:
        print("No, the conditions for a stabilizer code are not met.")
        global_answer = "No"
    
    return global_answer

if __name__ == '__main__':
    answer = check_stabilizer_code()
    # The final answer is wrapped in <<<>>>
    # print(f"\n<<<{answer}>>>") # This line is for the final output format, but we will just output the string directly.
