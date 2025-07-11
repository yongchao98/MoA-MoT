import sys

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage using dynamic programming.
    """
    # Problem data
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)
    num_knapsacks = len(capacities)

    # Sort items by weight to potentially improve pruning, though not strictly necessary
    # items = sorted(zip(weights, values))
    # weights = [item[0] for item in items]
    # values = [item[1] for item in items]


    # DP tables using dictionaries for sparse representation to save memory.
    # dp[(i, c0, c1, c2)] = max_value
    # decision[(i, c0, c1, c2)] = choice (0:skip, 1:K0, 2:K1, 3:K2)
    dp = {}
    decision = {}

    C0, C1, C2 = capacities[0], capacities[1], capacities[2]

    # Base case: Before considering any items, the value for any capacity combination is 0.
    for c0 in range(C0 + 1):
        for c1 in range(C1 + 1):
            for c2 in range(C2 + 1):
                dp[(0, c0, c1, c2)] = 0

    # Iterate through each item
    for i in range(1, num_items + 1):
        item_idx = i - 1
        v = values[item_idx]
        w = weights[item_idx]

        # Iterate through all possible capacity combinations
        for c0 in range(C0 + 1):
            for c1 in range(C1 + 1):
                for c2 in range(C2 + 1):
                    # Choice 0: Don't take the current item.
                    # The value is inherited from the state before considering this item.
                    best_val = dp.get((i - 1, c0, c1, c2), 0)
                    best_choice = 0

                    # Choice 1: Place in Knapsack 0
                    if c0 >= w:
                        val_k0 = v + dp.get((i - 1, c0 - w, c1, c2), 0)
                        if val_k0 > best_val:
                            best_val = val_k0
                            best_choice = 1

                    # Choice 2: Place in Knapsack 1
                    if c1 >= w:
                        val_k1 = v + dp.get((i - 1, c0, c1 - w, c2), 0)
                        if val_k1 > best_val:
                            best_val = val_k1
                            best_choice = 2

                    # Choice 3: Place in Knapsack 2
                    if c2 >= w:
                        val_k2 = v + dp.get((i - 1, c0, c1, c2 - w), 0)
                        if val_k2 > best_val:
                            best_val = val_k2
                            best_choice = 3

                    # Store the optimal value and decision for the current state
                    if best_val > 0:
                        dp[(i, c0, c1, c2)] = best_val
                        if best_choice != 0:
                            decision[(i, c0, c1, c2)] = best_choice
    
    # --- Backtracking to find the chosen items ---
    taken_item_values = []
    c0, c1, c2 = C0, C1, C2

    for i in range(num_items, 0, -1):
        choice = decision.get((i, c0, c1, c2), 0)
        item_idx = i - 1
        w = weights[item_idx]
        v = values[item_idx]

        if choice == 1: # Item was placed in Knapsack 0
            taken_item_values.append(v)
            c0 -= w
        elif choice == 2: # Item was placed in Knapsack 1
            taken_item_values.append(v)
            c1 -= w
        elif choice == 3: # Item was placed in Knapsack 2
            taken_item_values.append(v)
            c2 -= w
        # If choice is 0, item was not taken, capacities remain the same for the previous state.

    # --- Print the final result ---
    total_value = sum(taken_item_values)
    taken_item_values.sort(reverse=True)
    equation_str = " + ".join(map(str, taken_item_values))
    
    print(f"Maximum possible total value: {total_value}")
    print("Final Equation:")
    print(f"{total_value} = {equation_str}")


solve_multiple_knapsacks()
<<<569>>>