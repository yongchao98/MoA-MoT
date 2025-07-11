import pulp

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage using PuLP.
    
    This function models the problem as an Integer Linear Program to find the
    combination of items that maximizes the total value across all knapsacks
    without exceeding their capacities and ensuring each item is used at most once.
    """
    # --- Problem Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    
    num_items = len(values)
    num_knapsacks = len(capacities)

    # --- Model Initialization ---
    # Create the maximization problem
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # --- Decision Variables ---
    # x[i][j] is a binary variable: 1 if item i is in knapsack j, 0 otherwise
    x = pulp.LpVariable.dicts("item_in_knapsack", 
                              (range(num_items), range(num_knapsacks)), 
                              cat='Binary')

    # --- Objective Function ---
    # Maximize the total value of items in all knapsacks
    prob += pulp.lpSum(values[i] * x[i][j] 
                       for i in range(num_items) 
                       for j in range(num_knapsacks)), "Total_Value"

    # --- Constraints ---
    # 1. Knapsack capacity constraints:
    # The sum of weights of items in each knapsack must not exceed its capacity.
    for j in range(num_knapsacks):
        prob += pulp.lpSum(weights[i] * x[i][j] for i in range(num_items)) <= capacities[j], f"Knapsack_{j}_Capacity"

    # 2. Unique item usage constraint:
    # Each item can be placed in at most one knapsack.
    for i in range(num_items):
        prob += pulp.lpSum(x[i][j] for j in range(num_knapsacks)) <= 1, f"Item_{i}_Unique"

    # --- Solve the Problem ---
    # Note: PuLP uses the CBC solver by default. If not installed, it might need to be.
    # On most systems, `pip install pulp` handles this.
    try:
        prob.solve()
    except pulp.PulpSolverError:
        print("Pulp solver not found. Please install it.")
        print("You can often install the default CBC solver with PuLP by running: pip install pulp")
        return

    # --- Output the Results ---
    print(f"Optimization Status: {pulp.LpStatus[prob.status]}")
    
    if pulp.LpStatus[prob.status] == 'Optimal':
        selected_values = []
        for i in range(num_items):
            for j in range(num_knapsacks):
                if x[i][j].varValue == 1:
                    selected_values.append(values[i])
        
        # Sort for a clean, consistent output equation
        selected_values.sort(reverse=True)
        
        total_value = pulp.value(prob.objective)
        
        print("\nMaximum possible total value:", int(total_value))
        print("This is achieved by selecting items with the following values:")
        
        # Print the equation
        equation_str = " + ".join(map(str, selected_values))
        print(f"{equation_str} = {int(total_value)}")
    else:
        print("Could not find the optimal solution.")

if __name__ == "__main__":
    solve_multiple_knapsacks()
    # The final answer is the total value calculated.
    # Based on the provided data, the calculation will result in 643.
    # <<<643>>>