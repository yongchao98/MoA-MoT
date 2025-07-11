def solve_qubits():
    """
    Calculates the maximum number of logical qubits for a given surface code configuration.
    """
    # Number of separate surface code patches
    num_patches = 2

    # Number of holes in each patch
    num_holes_per_patch = 2

    # For a planar surface code, the number of logical qubits is equal to the number of holes.
    qubits_per_patch = num_holes_per_patch

    # The total number of logical qubits is the sum from all separate patches.
    total_logical_qubits = qubits_per_patch * num_patches

    print("Each patch can encode a number of logical qubits equal to its number of holes.")
    print(f"Number of logical qubits from Patch 1: {qubits_per_patch}")
    print(f"Number of logical qubits from Patch 2: {qubits_per_patch}")
    print("\nThe total number of logical qubits is the sum from both patches.")
    print(f"Final Equation: {qubits_per_patch} + {qubits_per_patch} = {total_logical_qubits}")

solve_qubits()