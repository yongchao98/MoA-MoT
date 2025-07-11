import numpy as np

def check_stabilizer_code():
    """
    Checks if a 4-qubit code can be described by the given stabilizers.
    Code: |0_L> = |0000>, |1_L> = |1111>
    Stabilizers: S1 = Z1*Z2, S2 = Z2*Z3, S3 = Z3*Z4
    """
    # --- Define basic states and operators ---
    # Pauli Z and Identity matrices
    I = np.identity(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Computational basis states |0> and |1>
    q0 = np.array([[1], [0]], dtype=complex)
    q1 = np.array([[0], [1]], dtype=complex)

    # --- Construct the 4-qubit logical states and stabilizers ---
    # Logical states |0_L> and |1_L>
    L0 = np.kron(np.kron(q0, q0), np.kron(q0, q0))
    L1 = np.kron(np.kron(q1, q1), np.kron(q1, q1))

    # Stabilizer operators
    S1 = np.kron(np.kron(Z, Z), np.kron(I, I))  # Z1 * Z2
    S2 = np.kron(np.kron(I, Z), np.kron(Z, I))  # Z2 * Z3
    S3 = np.kron(np.kron(I, I), np.kron(Z, Z))  # Z3 * Z4

    stabilizers = {'S1=Z1*Z2': S1, 'S2=Z2*Z3': S2, 'S3=Z3*Z4': S3}
    logical_states = {'|0_L>=|0000>': L0, '|1_L>=|1111>': L1}
    
    all_conditions_met = True

    print("--- Condition 1: Checking Commutativity of Stabilizers ---")
    
    # Check [S1, S2]
    comm_12 = S1 @ S2 - S2 @ S1
    if not np.allclose(comm_12, np.zeros((16, 16))):
        all_conditions_met = False
    print(f"[S1, S2] is zero: {np.allclose(comm_12, np.zeros((16, 16)))}")
    
    # Check [S1, S3]
    comm_13 = S1 @ S3 - S3 @ S1
    if not np.allclose(comm_13, np.zeros((16, 16))):
        all_conditions_met = False
    print(f"[S1, S3] is zero: {np.allclose(comm_13, np.zeros((16, 16)))}")

    # Check [S2, S3]
    comm_23 = S2 @ S3 - S3 @ S3
    if not np.allclose(comm_23, np.zeros((16, 16))):
        all_conditions_met = False
    print(f"[S2, S3] is zero: {np.allclose(comm_23, np.zeros((16, 16)))}\n")
    
    if not all_conditions_met:
        print("Stabilizers do not commute. This is not a valid stabilizer group.")
        return

    print("--- Condition 2: Checking Stabilization of Logical States ---")
    print("We check if S|psi> = +1 * |psi> for all stabilizers S and logical states |psi>.\n")

    # Check S1 = Z1*Z2
    print("Stabilizer S1 = Z1*Z2:")
    # On |0_L>
    print("  Applying to |0_L> = |0000>:")
    print("    S1|0000> = (Z1|0>)(Z2|0>)|00> = (+1|0>)(+1|0>)|00> = 1 * |0000>")
    res_0 = S1 @ L0
    is_stabilized_0 = np.allclose(res_0, L0)
    print(f"    Calculation correct: {is_stabilized_0}")
    if not is_stabilized_0: all_conditions_met = False
    # On |1_L>
    print("  Applying to |1_L> = |1111>:")
    print("    S1|1111> = (Z1|1>)(Z2|1>)|11> = (-1|1>)(-1|1>)|11> = 1 * |1111>")
    res_1 = S1 @ L1
    is_stabilized_1 = np.allclose(res_1, L1)
    print(f"    Calculation correct: {is_stabilized_1}\n")
    if not is_stabilized_1: all_conditions_met = False

    # Check S2 = Z2*Z3
    print("Stabilizer S2 = Z2*Z3:")
    # On |0_L>
    print("  Applying to |0_L> = |0000>:")
    print("    S2|0000> = |0>(Z2|0>)(Z3|0>)|0> = |0>(+1|0>)(+1|0>)|0> = 1 * |0000>")
    res_0 = S2 @ L0
    is_stabilized_0 = np.allclose(res_0, L0)
    print(f"    Calculation correct: {is_stabilized_0}")
    if not is_stabilized_0: all_conditions_met = False
    # On |1_L>
    print("  Applying to |1_L> = |1111>:")
    print("    S2|1111> = |1>(Z2|1>)(Z3|1>)|1> = |1>(-1|1>)(-1|1>)|1> = 1 * |1111>")
    res_1 = S2 @ L1
    is_stabilized_1 = np.allclose(res_1, L1)
    print(f"    Calculation correct: {is_stabilized_1}\n")
    if not is_stabilized_1: all_conditions_met = False

    # Check S3 = Z3*Z4
    print("Stabilizer S3 = Z3*Z4:")
    # On |0_L>
    print("  Applying to |0_L> = |0000>:")
    print("    S3|0000> = |00>(Z3|0>)(Z4|0>) = |00>(+1|0>)(+1|0>) = 1 * |0000>")
    res_0 = S3 @ L0
    is_stabilized_0 = np.allclose(res_0, L0)
    print(f"    Calculation correct: {is_stabilized_0}")
    if not is_stabilized_0: all_conditions_met = False
    # On |1_L>
    print("  Applying to |1_L> = |1111>:")
    print("    S3|1111> = |11>(Z3|1>)(Z4|1>) = |11>(-1|1>)(-1|1>) = 1 * |1111>")
    res_1 = S3 @ L1
    is_stabilized_1 = np.allclose(res_1, L1)
    print(f"    Calculation correct: {is_stabilized_1}\n")
    if not is_stabilized_1: all_conditions_met = False

    print("--- Conclusion ---")
    if all_conditions_met:
        print("Yes, all conditions are met.")
        print("The code can be considered a stabilizer code with the given stabilizers.")
        global final_answer
        final_answer = "Yes"
    else:
        print("No, one or more conditions failed.")
        print("The code cannot be considered a stabilizer code with the given stabilizers.")
        final_answer = "No"

if __name__ == '__main__':
    final_answer = ""
    check_stabilizer_code()
    print(f"<<<{final_answer}>>>")
