def calculate_logical_qubits():
    """
    Calculates the number of logical qubits for a surface code with multiple patches and holes.
    """
    # Number of separate surface code patches
    num_patches = 2

    # Number of holes in each patch
    num_holes_per_patch = 2

    # The formula for logical qubits (k) in a planar code with b boundaries (holes) is k = b - 1.
    # We calculate the number of logical qubits for a single patch first.
    qubits_per_patch = num_holes_per_patch - 1

    # The total number of logical qubits is the sum of qubits from each independent patch.
    total_qubits = num_patches * qubits_per_patch

    print(f"The number of logical qubits 'k' for a planar surface code with 'b' holes is given by the formula: k = b - 1.")
    print(f"For a single patch with {num_holes_per_patch} holes, the number of logical qubits is {num_holes_per_patch} - 1 = {qubits_per_patch}.")
    print(f"Since there are {num_patches} such patches, the total number of logical qubits is calculated as:")
    
    # Print the final equation with all numbers shown
    print(f"Final Equation: {num_patches} * ({num_holes_per_patch} - 1) = {total_qubits}")
    
    print(f"\nTherefore, at most {total_qubits} logical qubits can be encoded.")

calculate_logical_qubits()