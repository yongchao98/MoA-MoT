import numpy as np
from scipy.optimize import milp, LinearConstraint

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using Integer Linear Programming.
    """
    # 1. Define problem data
    values = np.array([24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80])
    weights = np.array([45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14])
    capacities = np.array([40, 120, 200])

    num_items = len(values)
    num_knapsacks = len(capacities)

    # 2. Define the ILP variables and objective function
    # Variables x_ij: 1 if item i is in knapsack j, 0 otherwise.
    # We flatten this into a 1D array of size num_items * num_knapsacks.
    # The variable for (item i, knapsack j) is at index i * num_knapsacks + j.
    # Objective: Maximize sum(values[i] * x_ij).
    # Since milp minimizes, we minimize -sum(values[i] * x_ij).
    # The coefficient for x_ij in the objective is -values[i].
    c = np.array([-v for v in values for _ in range(num_knapsacks)])

    # 3. Define the constraints in the form A @ x <= b
    # The matrix A has (num_knapsacks + num_items) rows and (num_items * num_knapsacks) columns.
    A_ub = np.zeros((num_knapsacks + num_items, num_items * num_knapsacks))
    b_ub = np.zeros(num_knapsacks + num_items)

    # Constraint Type 1: Knapsack capacity constraints
    # For each knapsack j: sum over i (weights[i] * x_ij) <= capacities[j]
    for j in range(num_knapsacks):
        for i in range(num_items):
            variable_index = i * num_knapsacks + j
            A_ub[j, variable_index] = weights[i]
        b_ub[j] = capacities[j]

    # Constraint Type 2: Unique item usage constraints
    # For each item i: sum over j (x_ij) <= 1
    for i in range(num_items):
        row_index = num_knapsacks + i
        for j in range(num_knapsacks):
            variable_index = i * num_knapsacks + j
            A_ub[row_index, variable_index] = 1
        b_ub[row_index] = 1

    # 4. Define integrality and bounds for variables
    # All variables are binary (integers between 0 and 1)
    integrality = np.ones_like(c)
    bounds = (0, 1)

    # Wrap constraints for SciPy's milp function
    constraints = LinearConstraint(A_ub, ub=b_ub)
    
    # 5. Solve the problem
    print("Solving the Multiple Knapsack Problem...")
    res = milp(c=c, constraints=constraints, integrality=integrality, bounds=bounds)

    # 6. Process and print the results
    if not res.success:
        print("Solver failed to find a solution.")
        print(f"Status: {res.status}, Message: {res.message}")
        return

    total_value = -res.fun
    assignments = np.round(res.x).reshape((num_items, num_knapsacks))

    print(f"\nMaximum total value: {total_value:.0f}")
    print("-" * 50)

    knapsack_values_list = [[] for _ in range(num_knapsacks)]
    knapsack_weights = [0] * num_knapsacks
    knapsack_values = [0] * num_knapsacks

    for i in range(num_items):
        for j in range(num_knapsacks):
            if assignments[i, j] > 0.5:
                knapsack_values_list[j].append(values[i])
                knapsack_weights[j] += weights[i]
                knapsack_values[j] += values[i]
    
    overall_value_sum_str_parts = []
    for j in range(num_knapsacks):
        # Sort values for consistent output
        knapsack_values_list[j].sort(reverse=True)
        value_list_str = " + ".join(map(str, knapsack_values_list[j]))
        if value_list_str: 
             overall_value_sum_str_parts.append(f"({value_list_str})")
        
        print(f"Knapsack {j+1} (Capacity: {capacities[j]})")
        print(f"  - Items (Values): {knapsack_values_list[j]}")
        print(f"  - Total Weight: {knapsack_weights[j]}")
        print(f"  - Total Value: {knapsack_values[j]}")
        if value_list_str:
             print(f"  - Equation: {knapsack_values[j]} = {value_list_str}")
        else:
             print(f"  - Equation: 0 = 0")
        print("-" * 50)
        
    print(f"Final Equation:")
    final_equation_str = " + ".join(overall_value_sum_str_parts)
    print(f"{total_value:.0f} = {final_equation_str}")


if __name__ == '__main__':
    solve_multiple_knapsack()
<<<733>>>