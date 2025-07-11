from ortools.sat.python import cp_model

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage using CP-SAT.
    """
    # --- Input Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # --- Create the Model ---
    model = cp_model.CpModel()

    # --- Create the Variables ---
    # x[i, j] is 1 if item i is put in knapsack j, and 0 otherwise.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # --- Define the Constraints ---
    # 1. Each item can be assigned to at most one knapsack.
    for i in all_items:
        model.AddAtMostOne(x[i, j] for j in all_knapsacks)

    # 2. The total weight in each knapsack cannot exceed its capacity.
    for j in all_knapsacks:
        model.Add(sum(weights[i] * x[i, j] for i in all_items) <= capacities[j])

    # --- Define the Objective Function ---
    # Maximize the total value of all packed items.
    objective_terms = []
    for i in all_items:
        for j in all_knapsacks:
            objective_terms.append(values[i] * x[i, j])
    model.Maximize(sum(objective_terms))

    # --- Solve the Model ---
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # --- Print the Results ---
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print(f"Maximum possible total value: {solver.ObjectiveValue()}\n")
        
        total_value_recalculated = 0
        chosen_values_for_equation = []

        for j in all_knapsacks:
            knapsack_weight = 0
            knapsack_value = 0
            print(f'Knapsack {j} (Capacity: {capacities[j]}):')
            for i in all_items:
                if solver.Value(x[i, j]) == 1:
                    print(f'  - Item {i} (Value: {values[i]}, Weight: {weights[i]})')
                    knapsack_weight += weights[i]
                    knapsack_value += values[i]
                    chosen_values_for_equation.append(values[i])
            print(f'  Packed Weight: {knapsack_weight}')
            print(f'  Packed Value: {knapsack_value}\n')
            total_value_recalculated += knapsack_value
        
        print("--- Final Equation ---")
        # Sort values for a consistent output order
        chosen_values_for_equation.sort(reverse=True)
        equation_str = ' + '.join(map(str, chosen_values_for_equation))
        print(f"{equation_str} = {total_value_recalculated}")

    else:
        print('No solution found.')

if __name__ == '__main__':
    solve_multiple_knapsacks()
<<<698>>>