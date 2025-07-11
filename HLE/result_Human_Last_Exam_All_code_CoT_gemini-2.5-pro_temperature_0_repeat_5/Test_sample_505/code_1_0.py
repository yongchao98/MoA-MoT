def solve_surface_code_qubits():
    """
    Calculates the maximum number of logical qubits for a given number of
    surface code patches, each with a specific number of holes.
    """
    # Number of separate surface code patches
    num_patches = 2

    # Number of holes in each patch
    num_holes_per_patch = 2

    # For a planar surface code, the number of logical qubits is equal to the number of holes.
    qubits_per_patch = num_holes_per_patch

    # The total number of logical qubits is the sum from all independent patches.
    total_logical_qubits = num_patches * qubits_per_patch

    print(f"Number of surface code patches: {num_patches}")
    print(f"Number of holes per patch: {num_holes_per_patch}")
    print(f"Logical qubits per patch = Number of holes = {qubits_per_patch}")
    print("\nSince there are two independent patches, the total number of logical qubits is the sum from each patch.")
    print(f"Total logical qubits = {qubits_per_patch} (from patch 1) + {qubits_per_patch} (from patch 2) = {total_logical_qubits}")

solve_surface_code_qubits()