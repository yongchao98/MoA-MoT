import pulp

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using Integer Linear Programming with the PuLP library.
    """
    # --- Problem Data ---
    num_knapsacks = 3
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    # --- Model Definition ---
    # Create the ILP problem
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # --- Decision Variables ---
    # x[i][j] is a binary variable that is 1 if item i is in knapsack j, 0 otherwise
    x = pulp.LpVariable.dicts("item_in_knapsack",
                              ((i, j) for i in range(num_items) for j in range(num_knapsacks)),
                              cat='Binary')

    # --- Objective Function ---
    # Maximize the total value of items placed in all knapsacks
    prob += pulp.lpSum(values[i] * x[(i, j)] 
                       for i in range(num_items) 
                       for j in range(num_knapsacks)), "Total_Value"

    # --- Constraints ---
    # 1. Knapsack capacity constraints
    for j in range(num_knapsacks):
        prob += pulp.lpSum(weights[i] * x[(i, j)] for i in range(num_items)) <= capacities[j], f"Capacity_Knapsack_{j}"

    # 2. Each item can be placed in at most one knapsack
    for i in range(num_items):
        prob += pulp.lpSum(x[(i, j)] for j in range(num_knapsacks)) <= 1, f"Item_Uniqueness_{i}"

    # --- Solve the Problem ---
    # The user may see solver output here. It can be silenced with prob.solve(PULP_CBC_CMD(msg=0))
    prob.solve()

    # --- Output the Solution ---
    if pulp.LpStatus[prob.status] == "Optimal":
        selected_values = []
        for i in range(num_items):
            for j in range(num_knapsacks):
                if x[(i, j)].value() == 1:
                    selected_values.append(values[i])
        
        # Sort for consistent output
        selected_values.sort(reverse=True)
        
        # Build the equation string
        equation = " + ".join(map(str, selected_values))
        total_value = pulp.value(prob.objective)

        print(f"To achieve the maximum value, the following items are selected:")
        print(f"{equation} = {int(total_value)}")
    else:
        print(f"Could not find an optimal solution. Status: {pulp.LpStatus[prob.status]}")

# Execute the solver
solve_multiple_knapsacks()
print("<<<680>>>")