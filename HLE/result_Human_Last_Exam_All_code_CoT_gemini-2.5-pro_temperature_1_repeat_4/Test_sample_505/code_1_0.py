def calculate_logical_qubits():
    """
    Calculates the maximum number of logical qubits for a given surface code configuration.
    """
    # Total number of holes available in the system.
    total_holes = 2
    # Total number of separate surface code patches.
    num_patches = 2

    print("The number of logical qubits (k) in a single patch of surface code is equal to the number of holes (h) it contains.")
    print("Formula per patch: k = h\n")

    # The total number of logical qubits is the sum of logical qubits from all patches,
    # which is equal to the total number of holes in the system.
    total_logical_qubits = total_holes

    print(f"With {num_patches} patches and a total of {total_holes} holes, the total number of logical qubits is {total_logical_qubits}.")
    print("This holds true regardless of how the holes are distributed. For example:\n")

    # Scenario 1: One hole in each of the two patches.
    patch1_holes_s1 = 1
    patch2_holes_s1 = 1
    total_qubits_s1 = patch1_holes_s1 + patch2_holes_s1
    print("Scenario 1 (1 hole per patch):")
    print(f"Logical Qubits = (Qubits in Patch 1) + (Qubits in Patch 2)")
    print(f"                 = {patch1_holes_s1} + {patch2_holes_s1} = {total_qubits_s1}\n")

    # Scenario 2: Both holes in one patch.
    patch1_holes_s2 = 2
    patch2_holes_s2 = 0
    total_qubits_s2 = patch1_holes_s2 + patch2_holes_s2
    print("Scenario 2 (2 holes in one patch):")
    print(f"Logical Qubits = (Qubits in Patch 1) + (Qubits in Patch 2)")
    print(f"                 = {patch1_holes_s2} + {patch2_holes_s2} = {total_qubits_s2}\n")

    print(f"Therefore, the maximum number of logical qubits that can be encoded is {total_logical_qubits}.")

calculate_logical_qubits()