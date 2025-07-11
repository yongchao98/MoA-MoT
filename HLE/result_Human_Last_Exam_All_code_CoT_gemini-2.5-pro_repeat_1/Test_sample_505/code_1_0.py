import math

def calculate_logical_qubits():
    """
    Calculates the total number of logical qubits for a given number of
    surface code patches, each with a specified number of holes.
    """
    # Number of independent surface code patches
    num_patches = 2

    # Number of holes created in each patch
    holes_per_patch = 2

    # For a single planar patch of surface code, the number of logical qubits
    # it can encode is equal to the number of holes.
    # Formula: k = h
    logical_qubits_per_patch = holes_per_patch

    # The total number of logical qubits is the sum from all independent patches.
    # Since they are identical, we multiply.
    total_logical_qubits = num_patches * logical_qubits_per_patch

    print("--- Calculation Breakdown ---")
    print(f"Number of surface code patches: {num_patches}")
    print(f"Number of holes per patch: {holes_per_patch}")
    print(f"Based on the formula k=h (logical qubits = holes), the logical qubits per patch is: {logical_qubits_per_patch}")
    print("\n--- Final Calculation ---")
    print(f"Total logical qubits = (Number of patches) * (Logical qubits per patch)")
    print(f"The final equation is: {num_patches} * {logical_qubits_per_patch} = {total_logical_qubits}")
    print(f"At most, {total_logical_qubits} logical qubits can be encoded.")

calculate_logical_qubits()