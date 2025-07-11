def solve_surface_code_qubits():
    """
    Calculates the maximum number of logical qubits for a given number
    of surface code patches and holes.
    """
    num_patches = 2
    total_holes = 2

    # The number of logical qubits (k) encoded in a single planar patch
    # is equal to the number of holes (h) in it. k = h.
    # For multiple patches, the total logical qubits is the sum of the
    # qubits from each patch.
    # k_total = h1 + h2 + ...
    #
    # Given h1 + h2 = total_holes, it follows that k_total = total_holes.
    # Therefore, the number of logical qubits is fixed by the total number of holes,
    # regardless of their distribution.

    max_logical_qubits = total_holes

    print(f"Configuration: {num_patches} patches of surface code and a total of {total_holes} holes.")
    print("Rule: The number of logical qubits encoded in a patch is equal to its number of holes.")
    print("\nLet's check all possible distributions of holes:")

    # Iterate through all possible ways to distribute the holes in the first patch
    for h1 in range(total_holes + 1):
        h2 = total_holes - h1
        
        # Calculate the total logical qubits for this distribution
        k1 = h1
        k2 = h2
        k_total = k1 + k2
        
        print(f"\n  Distribution:")
        print(f"    - Patch 1 has {h1} holes, encoding {k1} logical qubits.")
        print(f"    - Patch 2 has {h2} holes, encoding {k2} logical qubits.")
        print(f"  Equation for total logical qubits:")
        # Final equation with numbers, as requested
        print(f"    Total = {k1} + {k2} = {k_total}")

    print("\n" + "="*40)
    print(f"As shown, the total number of logical qubits is always {max_logical_qubits}.")
    print(f"Therefore, the maximum number of logical qubits is {max_logical_qubits}.")

solve_surface_code_qubits()