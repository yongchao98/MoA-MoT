def solve_reversal_moves():
    """
    Calculates the minimum moves to reverse a sequence of 100 elements
    with the given special swap operations.
    """
    
    # 1. Define problem parameters
    num_elements = 100
    # The non-adjacent swap is between elements at i and i+5.
    # This defines the size of the cyclic group for our classes.
    mod_class_size = 5
    
    # 2. Calculate the number of elements in each equivalence class
    elements_per_class = num_elements // mod_class_size
    
    print(f"The {num_elements} positions are partitioned into {mod_class_size} classes.")
    print(f"Each class contains {elements_per_class} elements that can be rearranged internally for free.\n")

    # 3. Determine the required permutation of classes for reversal
    # An element at position 'i' (1-indexed) is in class (i-1) % 5.
    # It must move to position '101-i'.
    # The new class is (101-i-1) % 5 = (100-i) % 5 = (-i) % 5.
    # Mapping k -> k':
    # if i = 5m+k+1, then k' = -(5m+k+1) % 5 = (-k-1) % 5 = (4-k) % 5.
    
    print("To reverse the sequence, the contents of the position classes must be permuted.")
    
    # This gives the permutation P = (0 4)(1 3)(2)
    # We calculate the cost for each cycle in the permutation.
    total_moves = 0
    
    # Cost for cycle (2): class 2 maps to itself
    cost_cycle_2 = 0
    print("Analyzing permutation cycle (2):")
    print(f"  - Class 2 contents must end up in class 2 positions.")
    print(f"  - Cost = {cost_cycle_2} moves.\n")
    total_moves += cost_cycle_2
    
    # Cost for cycle (0 4): swap contents of class 0 and 4
    # Distance between class 0 and 4 is min(|0-4|, 5 - |0-4|) = min(4, 1) = 1.
    # They are adjacent. Swapping them takes 'elements_per_class' moves.
    dist_0_4 = 1
    cost_cycle_0_4 = dist_0_4 * elements_per_class
    print("Analyzing permutation cycle (0 4):")
    print(f"  - Class 0 and 4 contents must be swapped.")
    print(f"  - These classes are adjacent (distance = {dist_0_4}).")
    print(f"  - Cost = {dist_0_4} * {elements_per_class} = {cost_cycle_0_4} moves.\n")
    total_moves += cost_cycle_0_4

    # Cost for cycle (1 3): swap contents of class 1 and 3
    # Distance between class 1 and 3 is min(|1-3|, 5 - |1-3|) = min(2, 3) = 2.
    # They are not adjacent. Swapping them using an intermediary class takes 3 steps.
    # Cost = 3 * elements_per_class
    dist_1_3 = 2
    cost_factor_dist_2 = 3 # To swap A,C via B: swap(A,B), swap(B,C), swap(A,B)
    cost_cycle_1_3 = cost_factor_dist_2 * elements_per_class
    print("Analyzing permutation cycle (1 3):")
    print(f"  - Class 1 and 3 contents must be swapped.")
    print(f"  - These classes are separated by 1 class (distance = {dist_1_3}).")
    print(f"  - This requires a 3-step swap using the intermediary class.")
    print(f"  - Cost = {cost_factor_dist_2} * {elements_per_class} = {cost_cycle_1_3} moves.\n")
    total_moves += cost_cycle_1_3
    
    # 4. Final Calculation
    print("The final equation for the total minimum number of moves is:")
    print(f"Total Moves = (Cost for cycle (2)) + (Cost for cycle (0 4)) + (Cost for cycle (1 3))")
    print(f"Total Moves = {cost_cycle_2} + {cost_cycle_0_4} + {cost_cycle_1_3}")
    print(f"Total Moves = {total_moves}")

solve_reversal_moves()