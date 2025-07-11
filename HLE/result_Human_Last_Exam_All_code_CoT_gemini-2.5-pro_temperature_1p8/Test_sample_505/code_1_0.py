def solve():
    """
    Calculates the maximum number of logical qubits that can be encoded in
    two patches of surface code, each with two holes.
    """

    # Number of independent patches
    num_patches = 2

    # Number of holes in each patch
    num_holes_per_patch = 2

    # For a single planar surface code patch with 'n' holes, the number of
    # logical qubits it can encode is equal to 'n'.
    # This is because the base patch provides 1 qubit and the 'n' holes
    # provide an additional 'n-1' qubits, for a total of 1 + (n-1) = n.
    qubits_per_patch = num_holes_per_patch

    # The total number of logical qubits is the sum of the qubits from each
    # independent patch. Since both patches are identical, we can multiply.
    total_qubits = num_patches * qubits_per_patch

    # Print the final calculation as an equation.
    # We construct the equation string to show the contribution from each patch.
    equation_parts = [str(qubits_per_patch) for _ in range(num_patches)]
    equation_str = " + ".join(equation_parts) + f" = {total_qubits}"
    
    print("The total number of logical qubits is the sum of the qubits from each patch.")
    print("Calculation:")
    print(equation_str)

solve()
<<<4>>>