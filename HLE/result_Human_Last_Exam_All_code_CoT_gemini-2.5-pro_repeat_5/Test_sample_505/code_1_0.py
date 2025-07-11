def solve_surface_code_qubits():
    """
    Calculates the maximum number of logical qubits that can be encoded in
    two patches of surface code, each with two holes.
    """
    # 1. Define the parameters based on the problem description.
    # We have two separate, identical systems (patches).
    num_patches = 2
    # A "patch" is a planar surface, so its genus (g) is 0.
    genus = 0
    # The problem states "two holes". To maximize the number of logical qubits,
    # we interpret this as each of the two patches having two holes (boundaries, b).
    num_holes_per_patch = 2

    # 2. Calculate the number of logical qubits for a single patch.
    # The formula is k = 2*g + b - 1.
    qubits_per_patch = 2 * genus + num_holes_per_patch - 1

    # 3. Calculate the total number of logical qubits for all patches.
    # Since the patches are independent, we sum the qubits from each one.
    total_logical_qubits = num_patches * qubits_per_patch

    # 4. Print the explanation and the final equation.
    print("To find the number of logical qubits, we use the formula k = 2g + b - 1 for each patch.")
    print(f"A single patch is a planar surface, so its genus (g) is {genus}.")
    print(f"Each patch has {num_holes_per_patch} holes (boundaries, b).")
    print("\nFirst, we calculate the logical qubits for one patch:")
    print(f"Qubits per patch = (2 * g) + b - 1")
    print(f"Qubits per patch = (2 * {genus}) + {num_holes_per_patch} - 1 = {qubits_per_patch}")

    print("\nSince there are two identical patches, the total is:")
    print(f"Total logical qubits = (Qubits per patch) * (Number of patches)")
    print(f"Total logical qubits = {qubits_per_patch} * {num_patches} = {total_logical_qubits}")

solve_surface_code_qubits()