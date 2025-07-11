def solve_qubits():
    """
    Calculates the maximum number of logical qubits that can be encoded in a
    system of separate surface code patches with a given number of holes.
    """
    # Number of separate surface code patches
    num_patches = 2

    # Total number of holes across all patches
    total_holes = 2

    # The formula for the total number of logical qubits (K) in a system
    # of P separate patches and H total holes is: K = H + P.
    # This is because each patch contributes 1 logical qubit, and each hole
    # adds an additional logical qubit.
    # K = (qubits from patches) + (qubits from holes)
    # K = P + H
    
    total_logical_qubits = total_holes + num_patches

    print("To find the maximum number of logical qubits, we use the formula:")
    print("Total Logical Qubits = (Total Number of Holes) + (Number of Patches)")
    print("\nGiven:")
    print(f"Number of Patches = {num_patches}")
    print(f"Total Number of Holes = {total_holes}")
    print("\nThe final calculation is:")
    print(f"{total_holes} + {num_patches} = {total_logical_qubits}")

solve_qubits()