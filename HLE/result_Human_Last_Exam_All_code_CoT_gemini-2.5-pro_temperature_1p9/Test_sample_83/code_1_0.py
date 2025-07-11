from ortools.sat.python import cp_model

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem using Google's OR-Tools CP-SAT solver.
    """
    # --- Input Data ---
    item_values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    item_weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    knapsack_capacities = [40, 120, 200]
    
    num_items = len(item_values)
    num_knapsacks = len(knapsack_capacities)

    # --- 1. Create the Model ---
    model = cp_model.CpModel()

    # --- 2. Create the Variables ---
    # x[i, j] is a boolean variable which is 1 if item j is placed in knapsack i.
    x = {}
    for i in range(num_knapsacks):
        for j in range(num_items):
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # --- 3. Define the Constraints ---
    # a) Each item can be assigned to at most one knapsack.
    for j in range(num_items):
        model.AddAtMostOne(x[i, j] for i in range(num_knapsacks))

    # b) The total weight in each knapsack cannot exceed its capacity.
    for i in range(num_knapsacks):
        model.Add(sum(x[i, j] * item_weights[j] for j in range(num_items)) <= knapsack_capacities[i])

    # --- 4. Define the Objective Function ---
    # Maximize the total value of all items placed in the knapsacks.
    objective_terms = []
    for i in range(num_knapsacks):
        for j in range(num_items):
            objective_terms.append(item_values[j] * x[i, j])
    model.Maximize(sum(objective_terms))

    # --- 5. Create a Solver and Solve the Model ---
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # --- 6. Print the Results ---
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        max_total_value = int(solver.ObjectiveValue())
        print(f"Maximum possible total value: {max_total_value}\n")
        
        selected_item_values = []
        # Display the contents of each knapsack and gather selected values.
        for i in range(num_knapsacks):
            knapsack_weight = 0
            knapsack_value = 0
            item_list = []
            for j in range(num_items):
                if solver.Value(x[i, j]) > 0.5:
                    item_list.append(f"  - Item (Value: {item_values[j]}, Weight: {item_weights[j]})")
                    knapsack_weight += item_weights[j]
                    knapsack_value += item_values[j]
                    selected_item_values.append(item_values[j])
            
            print(f"Knapsack {i+1} (Capacity: {knapsack_capacities[i]}, Used: {knapsack_weight}, Value: {knapsack_value}):")
            if not item_list:
                print("  - Empty")
            else:
                for item_str in item_list:
                    print(item_str)
            print("-" * 30)

        # Sort values for a clean equation
        selected_item_values.sort(reverse=True)
        equation = " + ".join(map(str, selected_item_values))
        print(f"\nFinal Equation:\n{equation} = {max_total_value}")

    else:
        print("The solver could not find an optimal solution.")

if __name__ == '__main__':
    solve_multiple_knapsack()

<<<692>>>