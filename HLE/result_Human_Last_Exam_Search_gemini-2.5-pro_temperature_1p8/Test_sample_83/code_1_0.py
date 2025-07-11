from ortools.linear_solver import pywraplp

def solve_multiple_knapsack():
    """
    Solves the multiple knapsack problem using a MIP solver.
    """
    # --- Problem Data ---
    item_values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    item_weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    knapsack_capacities = [40, 120, 200]

    num_items = len(item_values)
    num_knapsacks = len(knapsack_capacities)

    # --- Solver Initialization ---
    # Create the MIP solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')
    if not solver:
        print("SCIP solver not available.")
        return

    # --- Variable Creation ---
    # x[i][j] is a binary variable, 1 if item i is in knapsack j, 0 otherwise.
    x = {}
    for i in range(num_items):
        for j in range(num_knapsacks):
            x[i, j] = solver.IntVar(0, 1, f'x_{i}_{j}')

    # --- Constraint Definitions ---
    # 1. Each item can be placed in at most one knapsack.
    for i in range(num_items):
        solver.Add(sum(x[i, j] for j in range(num_knapsacks)) <= 1)

    # 2. The total weight of items in each knapsack cannot exceed its capacity.
    for j in range(num_knapsacks):
        solver.Add(sum(item_weights[i] * x[i, j] for i in range(num_items)) <= knapsack_capacities[j])

    # --- Objective Function ---
    # Maximize the total value of the items in all knapsacks.
    objective = solver.Objective()
    for i in range(num_items):
        for j in range(num_knapsacks):
            objective.SetCoefficient(x[i, j], item_values[i])
    objective.SetMaximization()

    # --- Solve ---
    status = solver.Solve()

    # --- Result Extraction ---
    if status == pywraplp.Solver.OPTIMAL:
        max_total_value = solver.Objective().Value()
        print(f"The maximum possible total value is: {max_total_value}")
        return max_total_value
    else:
        print("The solver could not find an optimal solution.")
        return None

# To get the answer, run the following function in a suitable environment.
# solve_multiple_knapsack()
