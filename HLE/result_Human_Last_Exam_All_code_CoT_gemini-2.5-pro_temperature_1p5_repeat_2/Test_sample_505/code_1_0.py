def solve_qubit_count():
    """
    Calculates the maximum number of logical qubits for the given surface code configuration.
    """
    # Number of holes in each of the two patches
    holes_per_patch = 2

    # For a single planar patch of surface code, the number of logical qubits
    # is equal to the number of holes.
    qubits_in_patch_1 = holes_per_patch
    qubits_in_patch_2 = holes_per_patch

    # The total number of logical qubits is the sum from the two separate patches.
    total_logical_qubits = qubits_in_patch_1 + qubits_in_patch_2

    print("The number of logical qubits in a single surface code patch is equal to the number of holes.")
    print(f"Number of qubits in the first patch = {qubits_in_patch_1}")
    print(f"Number of qubits in the second patch = {qubits_in_patch_2}")
    print("\nTotal logical qubits are the sum from both patches:")
    print(f"{qubits_in_patch_1} + {qubits_in_patch_2} = {total_logical_qubits}")

solve_qubit_count()