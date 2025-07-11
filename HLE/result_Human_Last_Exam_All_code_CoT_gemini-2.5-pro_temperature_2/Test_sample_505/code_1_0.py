def solve_qubit_count():
    """
    Calculates the maximum number of logical qubits for the given surface code configuration.
    """
    # Number of separate patches of surface code.
    num_patches = 2

    # The number of logical qubits in a single planar patch of surface code
    # is equal to the number of holes created in it.
    # The phrase "two patches of surface code with two holes" is ambiguous.
    # To find the maximum ("at most") number of logical qubits, we must assume
    # the configuration that maximizes this value. This occurs when each patch
    # has two holes.
    holes_per_patch = 2

    # Logical qubits per patch equals the number of holes in it.
    logical_qubits_per_patch = holes_per_patch

    # The total number of logical qubits for a system of multiple independent patches
    # is the sum of the logical qubits on each patch.
    total_logical_qubits = num_patches * logical_qubits_per_patch

    # Output the explanation and the final equation.
    print(f"The number of logical qubits on a planar surface code patch is equal to the number of holes, h.")
    print(f"The system has {num_patches} patches.")
    print(f"To find the maximum possible logical qubits, we assume each patch has {holes_per_patch} holes.")
    print(f"Total logical qubits = Number of patches * Holes per patch")
    print(f"{total_logical_qubits} = {num_patches} * {holes_per_patch}")

solve_qubit_count()