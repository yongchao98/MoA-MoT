def solve_qubits():
    """
    Calculates the maximum number of logical qubits that can be encoded in
    two patches of surface code, each with two holes.
    """
    # Number of separate surface code patches
    num_patches = 2

    # Number of holes (boundaries) in each patch
    num_holes_per_patch = 2

    # For a single planar surface code with 'b' boundaries (holes),
    # the number of logical qubits it can encode is k = b - 1.
    # We calculate the number of logical qubits for a single patch first.
    logical_qubits_per_patch = num_holes_per_patch - 1

    # The total number of logical qubits is the sum of qubits from each
    # independent patch. Since the patches are identical, we can multiply.
    total_logical_qubits = num_patches * logical_qubits_per_patch

    # Print the step-by-step calculation as requested
    print("Calculation for the number of logical qubits:")
    print(f"1. Logical qubits per patch = (Number of holes per patch) - 1")
    print(f"   Logical qubits per patch = {num_holes_per_patch} - 1 = {logical_qubits_per_patch}")
    print(f"2. Total logical qubits = (Number of patches) * (Logical qubits per patch)")
    print(f"   Total logical qubits = {num_patches} * {logical_qubits_per_patch} = {total_logical_qubits}")
    print(f"\nTherefore, two patches of surface code with two holes each can encode at most {total_logical_qubits} logical qubits.")

solve_qubits()