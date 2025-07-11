from ortools.sat.python import cp_model

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage.
    """
    # --- Problem Data ---
    item_values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    item_weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    knapsack_capacities = [40, 120, 200]

    num_items = len(item_values)
    num_knapsacks = len(knapsack_capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # --- Create the Model ---
    model = cp_model.CpModel()

    # --- Create the Decision Variables ---
    # x[i, j] is 1 if item i is put in knapsack j, and 0 otherwise.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = model.NewBoolVar(f'x_item_{i}_knapsack_{j}')

    # --- Define the Constraints ---
    # 1. Each item can be assigned to at most one knapsack.
    for i in all_items:
        model.AddAtMostOne(x[i, j] for j in all_knapsacks)

    # 2. The total weight in each knapsack cannot exceed its capacity.
    for j in all_knapsacks:
        model.Add(sum(item_weights[i] * x[i, j] for i in all_items) <= knapsack_capacities[j])

    # --- Define the Objective Function ---
    # Maximize the total value of all items in all knapsacks.
    objective_terms = []
    for i in all_items:
        for j in all_knapsacks:
            objective_terms.append(item_values[i] * x[i, j])
    model.Maximize(sum(objective_terms))

    # --- Solve the Model ---
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # --- Print the Results ---
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        max_value = solver.ObjectiveValue()
        print(f'Maximum possible total value: {max_value}\n')
        
        selected_values_for_equation = []

        for j in all_knapsacks:
            print("-" * 30)
            print(f'Knapsack {j+1} (Capacity: {knapsack_capacities[j]}):')
            knapsack_value = 0
            knapsack_weight = 0
            items_in_knapsack_details = []

            for i in all_items:
                if solver.Value(x[i, j]) == 1:
                    knapsack_value += item_values[i]
                    knapsack_weight += item_weights[i]
                    selected_values_for_equation.append(item_values[i])
                    items_in_knapsack_details.append(f"  - Item {i} (Value: {item_values[i]}, Weight: {item_weights[i]})")
            
            if items_in_knapsack_details:
                for detail in items_in_knapsack_details:
                    print(detail)
                print(f"\n  Total Value in Knapsack: {knapsack_value}")
                print(f"  Total Weight in Knapsack: {knapsack_weight} / {knapsack_capacities[j]}")
            else:
                print("  - Empty")
        
        print("-" * 30)
        
        # --- Print the final equation as requested ---
        print("\nFinal Equation for Total Value:")
        # Sort for a consistent, readable output
        selected_values_for_equation.sort(reverse=True)
        equation_str = " + ".join(map(str, selected_values_for_equation))
        print(f"{equation_str} = {int(max_value)}")

    else:
        print('No solution found for the given problem.')

if __name__ == '__main__':
    solve_multiple_knapsacks()
<<<650>>>