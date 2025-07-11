from ortools.sat.python import cp_model

def solve_multiple_knapsack():
    """Solves the multiple knapsack problem."""
    # Data
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)
    num_knapsacks = len(capacities)

    # Create the CP-SAT model.
    model = cp_model.CpModel()

    # Create the variables.
    # x[i, j] is a boolean variable, which is 1 if item i is placed in knapsack j.
    x = {}
    for i in range(num_items):
        for j in range(num_knapsacks):
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # Define the constraints.
    # Each item can be assigned to at most one knapsack.
    for i in range(num_items):
        model.AddAtMostOne(x[i, j] for j in range(num_knapsacks))

    # The total weight of the items in each knapsack must not exceed its capacity.
    for j in range(num_knapsacks):
        model.Add(sum(weights[i] * x[i, j] for i in range(num_items)) <= capacities[j])

    # Define the objective function.
    # Maximize the total value of the packed items.
    objective = []
    for i in range(num_items):
        for j in range(num_knapsacks):
            objective.append(values[i] * x[i, j])
    model.Maximize(sum(objective))

    # Create a solver and solve the model.
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # Print the solution.
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        total_value = solver.ObjectiveValue()
        # print(f'Total packed value: {total_value}')
        # total_weight = 0
        # for j in range(num_knapsacks):
        #     knapsack_weight = 0
        #     knapsack_value = 0
        #     print(f'Knapsack {j+1}:')
        #     for i in range(num_items):
        #         if solver.Value(x[i, j]) > 0:
        #             print(f'  Item {i} - Value: {values[i]}, Weight: {weights[i]}')
        #             knapsack_weight += weights[i]
        #             knapsack_value += values[i]
        #     print(f'  Total weight: {knapsack_weight} (Capacity: {capacities[j]})')
        #     print(f'  Total value: {knapsack_value}')
        #     total_weight += knapsack_weight
        return int(total_value)
    else:
        return "The problem does not have an optimal solution."

result = solve_multiple_knapsack()
# print(result)