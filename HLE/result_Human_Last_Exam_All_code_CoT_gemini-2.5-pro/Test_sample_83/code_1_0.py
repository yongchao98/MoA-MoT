from ortools.sat.python import cp_model

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage.
    This function finds the optimal combination of items to place in the
    knapsacks to maximize the total value, respecting capacity constraints
    and ensuring each item is used at most once.
    """
    # 1. Data Initialization
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    
    num_items = len(values)
    num_knapsacks = len(capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # 2. Create the CP-SAT Model
    model = cp_model.CpModel()

    # 3. Define Decision Variables
    # x[i, j] is a boolean variable that is true if item i is placed in knapsack j.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # 4. Define Constraints
    # Constraint 1: Each item can be placed in at most one knapsack.
    for i in all_items:
        model.AddAtMostOne(x[i, j] for j in all_knapsacks)

    # Constraint 2: The total weight of items in each knapsack cannot exceed its capacity.
    for j in all_knapsacks:
        knapsack_weight = sum(weights[i] * x[i, j] for i in all_items)
        model.Add(knapsack_weight <= capacities[j])

    # 5. Define the Objective Function
    # The objective is to maximize the total value of all packed items.
    total_value = sum(values[i] * x[i, j] for i in all_items for j in all_knapsacks)
    model.Maximize(total_value)

    # 6. Create a Solver and Solve the Model
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # 7. Extract and Print the Results
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        final_value = int(solver.ObjectiveValue())
        
        selected_item_values = []
        for i in all_items:
            for j in all_knapsacks:
                if solver.Value(x[i, j]) > 0:
                    selected_item_values.append(values[i])
        
        # Sort values for a consistent output format
        selected_item_values.sort(reverse=True)
        
        # Print the final equation
        equation_parts = [str(v) for v in selected_item_values]
        print(f"{' + '.join(equation_parts)} = {final_value}")

    else:
        print("No solution found for the problem.")

# Execute the function
solve_multiple_knapsacks()
<<<99 + 90 + 84 + 82 + 80 + 67 + 65 + 63 + 50 + 25 = 705>>>