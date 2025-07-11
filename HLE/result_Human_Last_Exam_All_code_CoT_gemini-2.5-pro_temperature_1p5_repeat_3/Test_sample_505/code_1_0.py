def calculate_logical_qubits():
    """
    Calculates the maximum number of logical qubits for a given number of
    surface code patches and holes.
    """
    # Number of separate patches of surface code
    num_patches = 2

    # Total number of holes across all patches
    num_holes = 2

    # Each patch contributes 1 logical qubit
    qubits_from_patches = num_patches * 1

    # Each hole adds the capacity for 1 additional logical qubit
    qubits_from_holes = num_holes * 1

    # The total is the sum of contributions from patches and holes
    total_logical_qubits = qubits_from_patches + qubits_from_holes

    # The final equation shows the contribution from patches and holes
    print(
        f"The maximum number of logical qubits is derived from {num_patches} patches and {num_holes} holes."
    )
    print(
        f"Total logical qubits = (Qubits from patches) + (Qubits from holes)"
    )
    print(
        f"Total logical qubits = {qubits_from_patches} + {qubits_from_holes} = {total_logical_qubits}"
    )


if __name__ == "__main__":
    calculate_logical_qubits()