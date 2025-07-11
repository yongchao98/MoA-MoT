def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using dynamic programming.
    """
    # --- Problem Definition ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    # --- Explanation of the Method ---
    print("This task is an instance of the Multiple Knapsack Problem (MKP).")
    print("The goal is to maximize the total value of items placed in three knapsacks, with the constraint that each item can be used at most once.")
    print("\nWe will use a dynamic programming approach to find the exact optimal solution.")
    print("The state dp[i][c1][c2][c3] will store the maximum value achievable using the first 'i' items with respective remaining knapsack capacities of c1, c2, and c3.")
    print("\nThe process can be computationally intensive but guarantees the best possible score.")
    print("\nCalculating the optimal value, please wait...\n")

    # --- Dynamic Programming Solution ---
    cap1, cap2, cap3 = capacities[0], capacities[1], capacities[2]
    
    # Initialize DP table
    # dp[i][w1][w2][w3] = max value using items 0..i-1 with capacities w1, w2, w3
    dp = [[[[0 for _ in range(cap3 + 1)] for _ in range(cap2 + 1)] for _ in range(cap1 + 1)] for _ in range(num_items + 1)]

    # Fill DP table
    for i in range(1, num_items + 1):
        item_index = i - 1
        v = values[item_index]
        w = weights[item_index]

        for c1 in range(cap1 + 1):
            for c2 in range(cap2 + 1):
                for c3 in range(cap3 + 1):
                    # Option 1: Don't take the item
                    value_no_place = dp[i-1][c1][c2][c3]
                    
                    options = [value_no_place]

                    # Option 2: Place in knapsack 1
                    if w <= c1:
                        options.append(dp[i-1][c1 - w][c2][c3] + v)
                    
                    # Option 3: Place in knapsack 2
                    if w <= c2:
                        options.append(dp[i-1][c1][c2 - w][c3] + v)
                    
                    # Option 4: Place in knapsack 3
                    if w <= c3:
                        options.append(dp[i-1][c1][c2][c3 - w] + v)
                    
                    dp[i][c1][c2][c3] = max(options)

    # The maximum value is the last cell of the DP table
    max_value = dp[num_items][cap1][cap2][cap3]

    # --- Reconstruct the Solution ---
    selected_values = []
    c1, c2, c3 = cap1, cap2, cap3

    for i in range(num_items, 0, -1):
        item_index = i - 1
        v = values[item_index]
        w = weights[item_index]

        # Check if the current value came from NOT placing the item
        if dp[i][c1][c2][c3] == dp[i-1][c1][c2][c3]:
            continue

        # Check if the item was placed in knapsack 1
        if w <= c1 and dp[i][c1][c2][c3] == dp[i-1][c1 - w][c2][c3] + v:
            selected_values.append(v)
            c1 -= w
            continue
            
        # Check if the item was placed in knapsack 2
        if w <= c2 and dp[i][c1][c2][c3] == dp[i-1][c1][c2 - w][c3] + v:
            selected_values.append(v)
            c2 -= w
            continue
            
        # Check if the item was placed in knapsack 3
        if w <= c3 and dp[i][c1][c2][c3] == dp[i-1][c1][c2][c3 - w] + v:
            selected_values.append(v)
            c3 -= w
            continue
            
    # --- Final Output ---
    selected_values.sort(reverse=True)
    equation = " + ".join(map(str, selected_values))
    
    print("Maximum possible total value found.")
    print("\nFinal Value Calculation:")
    print(f"{equation} = {max_value}")

# --- Execute the Solution ---
solve_multiple_knapsack()
<<<691>>>