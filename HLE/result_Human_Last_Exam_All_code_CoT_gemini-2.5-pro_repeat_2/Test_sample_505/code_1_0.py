def solve_surface_code_qubits():
    """
    Calculates the maximum number of logical qubits for a given surface code configuration.
    """
    # Parameters from the problem statement
    num_patches = 2
    total_holes = 2

    print(f"System Configuration:")
    print(f"- Number of separate surface code patches: {num_patches}")
    print(f"- Total number of holes across all patches: {total_holes}\n")

    print("Calculation Principle:")
    print("For a single planar surface code patch, the number of logical qubits it can encode is equal to its number of holes (k=h).\n")

    # The total number of logical qubits is the sum of qubits from each patch.
    # Total k = k_patch1 + k_patch2 = h_patch1 + h_patch2
    # Since h_patch1 + h_patch2 is the total number of holes, the total number of logical qubits
    # is simply equal to the total number of holes, regardless of distribution.
    total_logical_qubits = total_holes

    # To demonstrate, let's consider one possible distribution:
    # 1 hole in the first patch and 1 hole in the second patch.
    holes_in_patch1 = 1
    holes_in_patch2 = 1

    # Calculate logical qubits for each patch based on this distribution
    qubits_in_patch1 = holes_in_patch1
    qubits_in_patch2 = holes_in_patch2
    
    print("Illustrative Calculation (1 hole per patch):")
    print(f"Qubits in Patch 1 = {qubits_in_patch1}")
    print(f"Qubits in Patch 2 = {qubits_in_patch2}")
    
    # The final equation summing the parts
    print("\nFinal Equation:")
    print(f"Total Logical Qubits = {qubits_in_patch1} + {qubits_in_patch2} = {total_logical_qubits}")

solve_surface_code_qubits()