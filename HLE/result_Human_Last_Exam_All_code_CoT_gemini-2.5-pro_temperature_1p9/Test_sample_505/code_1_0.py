def solve_qubit_encoding():
    """
    Calculates the maximum number of logical qubits for the given surface code configuration.
    """
    # The number of separate patches of surface code.
    num_patches = 2

    # The number of holes in each individual patch.
    # In a planar surface code, the number of logical qubits equals the number of holes.
    holes_per_patch = 2
    logical_qubits_per_patch = holes_per_patch

    # The total number of logical qubits is the sum of the qubits from each independent patch.
    total_logical_qubits = num_patches * logical_qubits_per_patch

    # Print the explanation and the final equation.
    print(f"A single patch of planar surface code with {holes_per_patch} holes can encode {logical_qubits_per_patch} logical qubits.")
    print(f"Since there are {num_patches} independent patches, the total number of logical qubits is the sum from each patch.")
    print("\nFinal Equation:")
    print(f"{num_patches} (patches) * {logical_qubits_per_patch} (logical qubits per patch) = {total_logical_qubits} (total logical qubits)")

solve_qubit_encoding()