import numpy as np
from scipy.optimize import milp

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage using an ILP solver.
    """
    # --- Problem Data ---
    values = np.array([24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80])
    weights = np.array([45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14])
    capacities = np.array([40, 120, 200])
    n_items = len(values)
    n_knapsacks = len(capacities)

    # --- ILP Formulation ---
    # Decision variables: x_ij = 1 if item j is in knapsack i, 0 otherwise.
    # Total number of variables is n_knapsacks * n_items.
    # The variables are flattened into a single vector for the solver.
    n_vars = n_knapsacks * n_items

    # 1. Objective function: Maximize sum(v_j * x_ij)
    # The solver minimizes, so we minimize -sum(v_j * x_ij).
    # The objective coefficients vector is [-v_0, -v_1, ...,] repeated for each knapsack.
    c = -np.tile(values, n_knapsacks)

    # 2. Constraints (A_ub @ x <= b_ub)
    n_constraints = n_knapsacks + n_items
    A_ub = np.zeros((n_constraints, n_vars))
    b_ub = np.zeros(n_constraints)

    # Constraint type 1: Knapsack capacity constraints
    # For each knapsack i: sum(w_j * x_ij) <= C_i
    for i in range(n_knapsacks):
        start_col = i * n_items
        end_col = (i + 1) * n_items
        A_ub[i, start_col:end_col] = weights
        b_ub[i] = capacities[i]

    # Constraint type 2: Unique item usage constraints
    # For each item j: sum(x_ij over all knapsacks i) <= 1
    for j in range(n_items):
        constraint_row = n_knapsacks + j
        for i in range(n_knapsacks):
            # The variable for item j in knapsack i is at index (i * n_items + j)
            var_index = i * n_items + j
            A_ub[constraint_row, var_index] = 1
        b_ub[constraint_row] = 1

    # 3. Variable properties
    # All variables are binary (integers between 0 and 1).
    integrality = np.ones(n_vars)
    bounds = (0, 1)

    # --- Solve the ILP ---
    # Note: This requires scipy and numpy. You can install them with: pip install numpy scipy
    res = milp(c=c, A_ub=A_ub, b_ub=b_ub, integrality=integrality, bounds=bounds)

    # --- Process and Print Results ---
    if res.success:
        max_value = int(round(-res.fun))
        selected_vars = np.round(res.x).astype(int)

        selected_item_values = []
        # Check each item to see if it was placed in any knapsack
        for j in range(n_items):
            item_is_packed = sum(selected_vars[i * n_items + j] for i in range(n_knapsacks))
            if item_is_packed > 0:
                selected_item_values.append(values[j])

        # Sort values for a clean, consistent output
        selected_item_values.sort(reverse=True)

        # Print the final equation showing the sum of selected item values
        equation_parts = [str(v) for v in selected_item_values]
        print(f"{' + '.join(equation_parts)} = {max_value}")
    else:
        print("The solver could not find an optimal solution.")
        print(f"Solver status: {res.message}")

if __name__ == '__main__':
    solve_multiple_knapsacks()