import numpy as np
from scipy.optimize import milp, LinearConstraint

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using scipy's Mixed-Integer Linear Programming solver.
    """
    # --- Problem Data ---
    values = np.array([24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80])
    weights = np.array([45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14])
    capacities = np.array([40, 120, 200])
    n_knapsacks = len(capacities)
    n_items = len(values)

    # --- 1. Formulate the Objective Function ---
    # Decision variables x_ij: 1 if item j is in knapsack i, 0 otherwise.
    # We flatten these into a single vector of variables for the solver.
    # Objective: Maximize sum(value_j * x_ij).
    # Since milp minimizes, we minimize the negative of the sum.
    # The objective coefficients `c` are the item values, tiled for each knapsack, and negated.
    c = -np.tile(values, n_knapsacks)

    # --- 2. Formulate the Constraints ---
    # We have two sets of constraints:
    # A. Capacity: For each knapsack i, sum(weight_j * x_ij) <= capacity_i
    # B. Uniqueness: For each item j, sum(x_ij over all knapsacks i) <= 1
    
    num_vars = n_knapsacks * n_items
    num_constraints = n_knapsacks + n_items
    
    A = np.zeros((num_constraints, num_vars))
    b_u = np.zeros(num_constraints)

    # A. Capacity constraints
    for i in range(n_knapsacks):
        start_col = i * n_items
        end_col = (i + 1) * n_items
        A[i, start_col:end_col] = weights
        b_u[i] = capacities[i]

    # B. Uniqueness constraints
    for j in range(n_items):
        constraint_row = n_knapsacks + j
        for i in range(n_knapsacks):
            var_index = i * n_items + j
            A[constraint_row, var_index] = 1
        b_u[constraint_row] = 1

    # Combine into a single SciPy LinearConstraint object
    constraints = LinearConstraint(A, ub=b_u)

    # --- 3. Define Integrality and Bounds ---
    # All variables are binary (integer and between 0 and 1)
    integrality = np.ones_like(c)
    bounds = (0, 1)

    # --- 4. Solve the MILP Problem ---
    res = milp(c=c, constraints=constraints, integrality=integrality, bounds=bounds)

    # --- 5. Process and Print Results ---
    if not res.success:
        print("Optimization failed to find a solution.")
        print(f"Status: {res.status}, Message: {res.message}")
        return

    # Reshape the flat result vector into a (knapsack, item) matrix for easy interpretation
    assignments = np.round(res.x).reshape((n_knapsacks, n_items)).astype(int)

    print("Knapsack Assignments:\n")
    
    # Iterate through each knapsack to display its contents and calculate totals
    for i in range(n_knapsacks):
        # Find the indices of items assigned to the current knapsack
        assigned_item_indices = np.where(assignments[i, :] == 1)[0]
        
        knapsack_values = values[assigned_item_indices]
        knapsack_weights = weights[assigned_item_indices]
        
        knapsack_total_value = knapsack_values.sum()
        knapsack_total_weight = knapsack_weights.sum()
        
        # Format the equation string showing the sum of values
        value_str = " + ".join(map(str, knapsack_values))
        if not value_str:
            value_str = "0"
            
        print(f"Knapsack {i+1} (Capacity: {capacities[i]}, Total Weight: {knapsack_total_weight}):")
        print(f"{value_str} = {knapsack_total_value}\n")

    # The optimal value is the negative of the minimized objective function result
    total_value = int(round(-res.fun))
    print("----------------------------------------")
    print(f"Maximum Total Value: {total_value}")

    print(f"\n<<<668>>>")

if __name__ == '__main__':
    solve_multiple_knapsacks()