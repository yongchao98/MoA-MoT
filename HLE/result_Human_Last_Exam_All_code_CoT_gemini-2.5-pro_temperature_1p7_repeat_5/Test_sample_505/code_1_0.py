def solve_qubit_count():
    """
    Calculates the maximum number of logical qubits in a system of
    two surface code patches, each with two holes.
    """

    # Number of independent surface code patches
    num_patches = 2

    # Number of holes in each patch
    num_holes_per_patch = 2

    # For a planar surface code, the number of logical qubits (k) that can be
    # encoded in a single patch is equal to the number of holes (h).
    logical_qubits_per_patch = num_holes_per_patch

    # The total number of logical qubits is the sum from all independent patches.
    total_logical_qubits = num_patches * logical_qubits_per_patch

    print("To find the total number of logical qubits, we first find the number per patch and then sum them.")
    print(f"Number of logical qubits per patch = Number of holes = {logical_qubits_per_patch}")
    print("Total logical qubits = (Qubits in Patch 1) + (Qubits in Patch 2)")
    print(f"Total logical qubits = {logical_qubits_per_patch} + {logical_qubits_per_patch} = {total_logical_qubits}")

solve_qubit_count()