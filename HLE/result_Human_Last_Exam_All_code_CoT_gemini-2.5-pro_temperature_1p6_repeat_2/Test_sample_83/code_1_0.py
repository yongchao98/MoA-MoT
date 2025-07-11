import sys
try:
    from ortools.sat.python import cp_model
except ImportError:
    print("ortools library not found. Please install it with 'pip install ortools'")
    sys.exit(1)

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage using ortools.
    """
    # --- Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(capacities)

    # --- Model Creation ---
    model = cp_model.CpModel()

    # --- Decision Variables ---
    # x[i][j] is 1 if item i is put in knapsack j, and 0 otherwise.
    x = {}
    for i in range(num_items):
        for j in range(num_knapsacks):
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # --- Constraints ---
    # 1. Each item can be assigned to at most one knapsack.
    for i in range(num_items):
        model.AddAtMostOne(x[i, j] for j in range(num_knapsacks))

    # 2. The total weight in each knapsack cannot exceed its capacity.
    for j in range(num_knapsacks):
        model.Add(sum(weights[i] * x[i, j] for i in range(num_items)) <= capacities[j])

    # --- Objective Function ---
    # Maximize the total value of all items in all knapsacks.
    objective_terms = []
    for i in range(num_items):
        for j in range(num_knapsacks):
            objective_terms.append(values[i] * x[i, j])
    model.Maximize(sum(objective_terms))

    # --- Solver ---
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # --- Output Results ---
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print(f'Maximum possible total value: {solver.ObjectiveValue():.0f}\n')
        total_weight = 0
        total_value = 0
        selected_item_values = []

        for j in range(num_knapsacks):
            knapsack_weight = 0
            knapsack_value = 0
            packed_items = []
            print(f'--- Knapsack {j+1} (Capacity: {capacities[j]}) ---')
            for i in range(num_items):
                if solver.Value(x[i, j]) == 1:
                    packed_items.append(f"Item {i+1} (Value: {values[i]}, Weight: {weights[i]})")
                    knapsack_weight += weights[i]
                    knapsack_value += values[i]
                    selected_item_values.append(values[i])
            if not packed_items:
                print("  (empty)")
            else:
                for item_str in packed_items:
                  print(f"  {item_str}")
            print(f"  Total weight: {knapsack_weight}")
            print(f"  Total value: {knapsack_value}\n")
            total_weight += knapsack_weight
            total_value += knapsack_value
        
        # Sort values for consistent output equation
        selected_item_values.sort(reverse=True)
        
        # Print the final equation
        equation_str = " + ".join(map(str, selected_item_values))
        print("Final value equation:")
        print(f"{equation_str} = {total_value}")

        print(f"\n<<<668>>>")
    else:
        print('No solution found.')

if __name__ == '__main__':
    solve_multiple_knapsacks()