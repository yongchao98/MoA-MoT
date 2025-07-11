import numpy as np
from scipy.optimize import milp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using an Integer Linear Programming approach.
    """
    # --- Problem Data ---
    values = np.array([24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80])
    weights = np.array([45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14])
    capacities = np.array([40, 120, 200])
    num_knapsacks = len(capacities)
    num_items = len(values)

    # --- ILP Model Formulation ---
    # Decision variables x_ij = 1 if item i is in knapsack j, 0 otherwise.
    # Variables are flattened into a 1D array: [x_00, x_01, x_02, x_10, x_11, x_12, ...]
    num_vars = num_items * num_knapsacks

    # Objective function: Maximize sum(v_i * x_ij)
    # The solver minimizes, so we use negative values: Minimize sum(-v_i * x_ij)
    c = -np.array([v for v in values for _ in range(num_knapsacks)])

    # Constraint matrix A_ub and vector b_ub for A_ub @ x <= b_ub
    A_ub = np.zeros((num_knapsacks + num_items, num_vars))
    b_ub = np.zeros(num_knapsacks + num_items)

    # 1. Knapsack capacity constraints: sum(w_i * x_ij) <= C_j for each knapsack j
    for j in range(num_knapsacks):
        b_ub[j] = capacities[j]
        for i in range(num_items):
            # The coefficient for variable x_ij is the weight of item i
            A_ub[j, i * num_knapsacks + j] = weights[i]

    # 2. Item uniqueness constraints: sum(x_ij over j) <= 1 for each item i
    for i in range(num_items):
        row_idx = num_knapsacks + i
        b_ub[row_idx] = 1
        for j in range(num_knapsacks):
            # The coefficients for variables x_i0, x_i1, ... are 1
            A_ub[row_idx, i * num_knapsacks + j] = 1

    # All variables are binary (Integer, with bounds 0 and 1)
    integrality = np.ones(num_vars)
    bounds = (0, 1) # Scipy broadcasts this to all variables

    # --- Solve the ILP ---
    print("Solving the ILP problem... This may take a moment.")
    result = milp(c=c, A_ub=A_ub, b_ub=b_ub, integrality=integrality, bounds=bounds)

    # --- Process and Print Results ---
    if result.success:
        max_value = -result.fun
        assignments = np.round(result.x).reshape((num_items, num_knapsacks))
        
        selected_item_values = []
        for i in range(num_items):
            # Check if item i was placed in any knapsack
            if assignments[i, :].sum() > 0:
                selected_item_values.append(values[i])
                
        # Sort values for a consistent, readable equation
        selected_item_values.sort(reverse=True)
        
        total_value_from_items = sum(selected_item_values)

        print(f"\nMaximum possible total value: {int(round(max_value))}")
        print("\nFinal Equation:")
        
        # Output the equation with each number as requested
        equation_str = " + ".join(map(str, selected_item_values))
        print(f"{equation_str} = {int(round(total_value_from_items))}")

    else:
        print("\nCould not find an optimal solution.")
        print(f"Solver status: {result.status} - {result.message}")

if __name__ == '__main__':
    solve_multiple_knapsack()