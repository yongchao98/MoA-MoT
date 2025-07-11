import numpy as np

def check_stabilizer_code():
    """
    Checks if the 4-qubit code with |0_L> = |0000> and |1_L> = |1111>
    is a stabilizer code for the stabilizers Z1Z2, Z2Z3, and Z3Z4.
    """
    # Define Pauli Z and Identity matrices for computational verification,
    # though the explanation will be based on the analytical results.
    I = np.array([[1, 0], [0, 1]])
    Z = np.array([[1, 0], [0, -1]])
    
    # Define logical states as vectors
    q0 = np.array([1, 0])
    q1 = np.array([0, 1])
    L0 = np.kron(np.kron(np.kron(q0, q0), q0), q0)
    L1 = np.kron(np.kron(np.kron(q1, q1), q1), q1)
    
    # Define stabilizer generators as matrices
    S1 = np.kron(np.kron(np.kron(Z, Z), I), I)
    S2 = np.kron(np.kron(np.kron(I, Z), Z), I)
    S3 = np.kron(np.kron(np.kron(I, I), Z), Z)

    print("To be a stabilizer code, all logical states must be eigenvectors with eigenvalue +1")
    print("for all stabilizer generators. This means S|ψ> = 1 * |ψ>.")
    print("We will check this for both |0_L> and |1_L>.\n")
    print("Recall: Z|0> = +1 * |0>  and  Z|1> = -1 * |1>\n")
    print("-" * 75)

    # --- Check |0_L> ---
    print("Checking stabilizer actions on |0_L> = |0000>:")

    # S1 on |0_L>
    print("\n1. Action of S1 = Z1*Z2:")
    print("   S1 |0_L> = (Z1*Z2) |0000> = (Z|0>)_1 (Z|0>)_2 |0>_3 |0>_4")
    print("   S1 |0_L> = (+1)|0>_1 * (+1)|0>_2 * |0>_3 * |0>_4 = 1 * |0000> = |0_L>")
    is_stabilized_s1_l0 = np.allclose(S1 @ L0, L0)
    
    # S2 on |0_L>
    print("\n2. Action of S2 = Z2*Z3:")
    print("   S2 |0_L> = (Z2*Z3) |0000> = |0>_1 (Z|0>)_2 (Z|0>)_3 |0>_4")
    print("   S2 |0_L> = |0>_1 * (+1)|0>_2 * (+1)|0>_3 * |0>_4 = 1 * |0000> = |0_L>")
    is_stabilized_s2_l0 = np.allclose(S2 @ L0, L0)

    # S3 on |0_L>
    print("\n3. Action of S3 = Z3*Z4:")
    print("   S3 |0_L> = (Z3*Z4) |0000> = |0>_1 |0>_2 (Z|0>)_3 (Z|0>)_4")
    print("   S3 |0_L> = |0>_1 * |0>_2 * (+1)|0>_3 * (+1)|0>_4 = 1 * |0000> = |0_L>")
    is_stabilized_s3_l0 = np.allclose(S3 @ L0, L0)

    print("\n=> Result: |0_L> is stabilized by all three generators.")
    print("-" * 75)

    # --- Check |1_L> ---
    print("Checking stabilizer actions on |1_L> = |1111>:")

    # S1 on |1_L>
    print("\n1. Action of S1 = Z1*Z2:")
    print("   S1 |1_L> = (Z1*Z2) |1111> = (Z|1>)_1 (Z|1>)_2 |1>_3 |1>_4")
    print("   S1 |1_L> = (-1)|1>_1 * (-1)|1>_2 * |1>_3 * |1>_4 = 1 * |1111> = |1_L>")
    is_stabilized_s1_l1 = np.allclose(S1 @ L1, L1)

    # S2 on |1_L>
    print("\n2. Action of S2 = Z2*Z3:")
    print("   S2 |1_L> = (Z2*Z3) |1111> = |1>_1 (Z|1>)_2 (Z|1>)_3 |1>_4")
    print("   S2 |1_L> = |1>_1 * (-1)|1>_2 * (-1)|1>_3 * |1>_4 = 1 * |1111> = |1_L>")
    is_stabilized_s2_l1 = np.allclose(S2 @ L1, L1)

    # S3 on |1_L>
    print("\n3. Action of S3 = Z3*Z4:")
    print("   S3 |1_L> = (Z3*Z4) |1111> = |1>_1 |1>_2 (Z|1>)_3 (Z|1>)_4")
    print("   S3 |1_L> = |1>_1 * |1>_2 * (-1)|1>_3 * (-1)|1>_4 = 1 * |1111> = |1_L>")
    is_stabilized_s3_l1 = np.allclose(S3 @ L1, L1)
    
    print("\n=> Result: |1_L> is stabilized by all three generators.")
    print("-" * 75)

    # Final Conclusion
    all_stabilized = (is_stabilized_s1_l0 and is_stabilized_s2_l0 and is_stabilized_s3_l0 and
                      is_stabilized_s1_l1 and is_stabilized_s2_l1 and is_stabilized_s3_l1)

    if all_stabilized:
        print("\nFinal Conclusion: YES.")
        print("Since both logical basis states, |0_L> and |1_L>, are eigenvectors with eigenvalue +1")
        print("for all proposed generators, the code can be considered a stabilizer code with")
        print("stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.")
    else:
        print("\nFinal Conclusion: NO.")
        print("At least one logical basis state was not stabilized by one of the generators.")

if __name__ == "__main__":
    check_stabilizer_code()
