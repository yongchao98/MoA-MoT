def check_stabilizer_code():
    """
    Checks if a 4-qubit code is a stabilizer code for the given stabilizers.
    The code prints the step-by-step verification.
    """
    print("We need to check if the logical states are stabilized by the operators,")
    print("meaning S|psi_L> = +1 * |psi_L> for all stabilizers S and logical states |psi_L>.")
    print("\nFirst, let's recall the action of the Pauli Z operator:")
    print("Z|0> = 1 * |0>")
    print("Z|1> = -1 * |1>")

    print("\n--- Step 1: Checking Stabilizer S1 = Z1*Z2 ---")
    print("\nOn state |0_L> = |0000>:")
    print("  S1 |0_L> = (Z1*Z2) |0000> = (Z1|0>) (Z2|0>) |0> |0>")
    print("           = (1*|0>) (1*|0>) |0> |0>")
    print("           = 1 * |0000>")
    print("Result: |0_L> is stabilized by S1.")

    print("\nOn state |1_L> = |1111>:")
    print("  S1 |1_L> = (Z1*Z2) |1111> = (Z1|1>) (Z2|1>) |1> |1>")
    print("           = ((-1)*|1>) ((-1)*|1>) |1> |1>")
    print("           = (-1)*(-1) * |1111>")
    print("           = 1 * |1111>")
    print("Result: |1_L> is stabilized by S1.")

    print("\n--- Step 2: Checking Stabilizer S2 = Z2*Z3 ---")
    print("\nOn state |0_L> = |0000>:")
    print("  S2 |0_L> = (Z2*Z3) |0000> = |0> (Z2|0>) (Z3|0>) |0>")
    print("           = |0> (1*|0>) (1*|0>) |0>")
    print("           = 1 * |0000>")
    print("Result: |0_L> is stabilized by S2.")

    print("\nOn state |1_L> = |1111>:")
    print("  S2 |1_L> = (Z2*Z3) |1111> = |1> (Z2|1>) (Z3|1>) |1>")
    print("           = |1> ((-1)*|1>) ((-1)*|1>) |1>")
    print("           = (-1)*(-1) * |1111>")
    print("           = 1 * |1111>")
    print("Result: |1_L> is stabilized by S2.")

    print("\n--- Step 3: Checking Stabilizer S3 = Z3*Z4 ---")
    print("\nOn state |0_L> = |0000>:")
    print("  S3 |0_L> = (Z3*Z4) |0000> = |0> |0> (Z3|0>) (Z4|0>)")
    print("           = |0> |0> (1*|0>) (1*|0>)")
    print("           = 1 * |0000>")
    print("Result: |0_L> is stabilized by S3.")

    print("\nOn state |1_L> = |1111>:")
    print("  S3 |1_L> = (Z3*Z4) |1111> = |1> |1> (Z3|1>) (Z4|1>)")
    print("           = |1> |1> ((-1)*|1>) ((-1)*|1>)")
    print("           = (-1)*(-1) * |1111>")
    print("           = 1 * |1111>")
    print("Result: |1_L> is stabilized by S3.")

    print("\n--- Final Conclusion ---")
    print("Since both logical basis states, |0_L> and |1_L>, are eigenvectors with eigenvalue +1 for all three operators,")
    print("the code can indeed be considered a stabilizer code with stabilizers Z1*Z2, Z2*Z3, and Z3*Z4.")

if __name__ == '__main__':
    check_stabilizer_code()