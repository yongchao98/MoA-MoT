def solve_logical_qubits():
    """
    Calculates the maximum number of logical qubits for the given surface code setup.
    """
    # 1. Define the parameters based on the problem statement.
    num_patches = 2
    # "two holes" is interpreted as two inner boundaries (defects) per patch.
    num_inner_holes_per_patch = 2

    # 2. Calculate the total number of boundaries for a single patch.
    # Each patch has 1 outer boundary plus the specified inner holes.
    num_boundaries_per_patch = num_inner_holes_per_patch + 1

    # 3. Calculate the number of logical qubits for a single patch.
    # The formula for a planar code is k = b - 1, where b is the number of boundaries.
    logical_qubits_per_patch = num_boundaries_per_patch - 1

    # 4. Calculate the total number of logical qubits for all patches.
    # Since the patches are independent, we sum the logical qubits from each.
    total_logical_qubits = num_patches * logical_qubits_per_patch

    # 5. Print the step-by-step reasoning and the final result.
    print("To find the maximum number of logical qubits, we follow these steps:")
    print(f"1. We have {num_patches} separate patches of surface code.")
    print(f"2. Each patch is described as having {num_inner_holes_per_patch} holes. In a planar code, a patch inherently has 1 outer boundary. The 'holes' are additional inner boundaries.")
    print(f"3. Thus, each patch has 1 (outer boundary) + {num_inner_holes_per_patch} (inner holes) = {num_boundaries_per_patch} boundaries in total.")
    print(f"4. The number of logical qubits in a single planar patch is given by the formula: (Number of Boundaries) - 1.")
    print(f"   So, for one patch: {num_boundaries_per_patch} - 1 = {logical_qubits_per_patch} logical qubits.")
    print(f"5. Since the two patches are independent, the total number of logical qubits is the sum from both.")
    print("\nFinal calculation:")
    print(f"Total Logical Qubits = (Qubits per Patch) * (Number of Patches)")
    print(f"Total Logical Qubits = {logical_qubits_per_patch} * {num_patches} = {total_logical_qubits}")

solve_logical_qubits()