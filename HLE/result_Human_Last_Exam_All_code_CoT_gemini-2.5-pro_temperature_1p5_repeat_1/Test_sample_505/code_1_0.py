def solve_logical_qubits():
    """
    Calculates the maximum number of logical qubits for a given number
    of surface code patches and holes.
    """
    num_patches = 2
    total_holes = 2

    print("The number of logical qubits in a single planar surface code patch is equal to the number of holes it contains.")
    print("The total number of logical qubits is the sum of the qubits from all individual patches.")
    print(f"\nWe have {num_patches} patches and a total of {total_holes} holes.")
    print("Let's consider the possible distributions of holes:")

    # Case 1: 2 holes in the first patch, 0 in the second.
    holes_patch1_case1 = 2
    holes_patch2_case1 = 0
    logical_qubits_case1 = holes_patch1_case1 + holes_patch2_case1
    
    print("\nDistribution 1: One patch has 2 holes, the other has 0.")
    print(f"Total Logical Qubits = (qubits from patch 1) + (qubits from patch 2)")
    print(f"Total Logical Qubits = {holes_patch1_case1} + {holes_patch2_case1} = {logical_qubits_case1}")

    # Case 2: 1 hole in the first patch, 1 in the second.
    holes_patch1_case2 = 1
    holes_patch2_case2 = 1
    logical_qubits_case2 = holes_patch1_case2 + holes_patch2_case2

    print("\nDistribution 2: Each patch has 1 hole.")
    print(f"Total Logical Qubits = (qubits from patch 1) + (qubits from patch 2)")
    print(f"Total Logical Qubits = {holes_patch1_case2} + {holes_patch2_case2} = {logical_qubits_case2}")

    # The maximum is the same in all cases.
    max_logical_qubits = total_holes
    print(f"\nSince the total is always {max_logical_qubits}, this is the maximum possible value.")

solve_logical_qubits()