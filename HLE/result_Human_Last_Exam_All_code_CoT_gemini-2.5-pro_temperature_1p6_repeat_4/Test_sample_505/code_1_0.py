def solve_logical_qubits():
    """
    Calculates the number of logical qubits for a given surface code configuration.
    """
    # Define the parameters of the problem
    num_patches = 2
    num_holes_per_patch = 2
    
    # For a planar surface code patch (genus g=0), the number of logical qubits (k)
    # is equal to the number of holes (h).
    # This comes from the general formula k = 2g + b - 1,
    # where g=0 (planar) and b = h + 1 (number of boundaries).
    # So, k = 2*0 + (h + 1) - 1 = h.
    k_per_patch = num_holes_per_patch
    
    # The total number of logical qubits is the sum from all independent patches.
    total_k = num_patches * k_per_patch
    
    print(f"Configuration:")
    print(f"Number of separate surface code patches: {num_patches}")
    print(f"Number of holes in each patch: {num_holes_per_patch}")
    print("-" * 30)
    print(f"Calculation per patch:")
    print(f"Logical qubits per patch = Number of holes per patch")
    print(f"Logical qubits per patch = {k_per_patch}")
    print("-" * 30)
    print(f"Total Calculation:")
    print(f"Total logical qubits = (Number of patches) * (Logical qubits per patch)")
    # Final equation as requested
    print(f"Total logical qubits = {num_patches} * {k_per_patch} = {total_k}")

solve_logical_qubits()