# Note: This code requires the 'ortools' library.
# You can install it by running: pip install ortools

from ortools.sat.python import cp_model

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage using a CP-SAT solver.
    """
    # 1. Define Problem Data
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # 2. Create the CP-SAT Model
    model = cp_model.CpModel()

    # 3. Define the Decision Variables
    # x[i, j] is a boolean variable, true if item i is placed in knapsack j.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # 4. Define the Constraints
    # a) Each item can be assigned to at most one knapsack.
    for i in all_items:
        model.AddAtMostOne(x[i, j] for j in all_knapsacks)

    # b) The total weight in each knapsack must not exceed its capacity.
    for j in all_knapsacks:
        model.Add(sum(x[i, j] * weights[i] for i in all_items) <= capacities[j])

    # 5. Define the Objective Function
    # Maximize the total value of all packed items.
    objective_terms = []
    for i in all_items:
        for j in all_knapsacks:
            objective_terms.append(x[i, j] * values[i])
    model.Maximize(sum(objective_terms))

    # 6. Create the Solver and Solve the Model
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # 7. Print the Solution
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        total_value = int(solver.ObjectiveValue())
        print(f'Maximum possible total value: {total_value}\n')

        packed_item_values = []
        for j in all_knapsacks:
            knapsack_value = 0
            knapsack_weight = 0
            print(f'Knapsack {j + 1} (Capacity: {capacities[j]})')
            
            items_in_knapsack = []
            for i in all_items:
                if solver.Value(x[i, j]) > 0:
                    items_in_knapsack.append((values[i], weights[i], i + 1))
                    packed_item_values.append(values[i])
                    knapsack_value += values[i]
                    knapsack_weight += weights[i]
            
            # Sort items by value for cleaner display
            items_in_knapsack.sort(key=lambda item: item[0], reverse=True)
            for value, weight, item_num in items_in_knapsack:
                 print(f'  - Item {item_num} (Value: {value}, Weight: {weight})')

            print(f'Total value for this knapsack: {knapsack_value}')
            print(f'Total weight for this knapsack: {knapsack_weight}\n')

        # Sort values for a consistent final equation
        packed_item_values.sort(reverse=True)
        equation_str = " + ".join(map(str, packed_item_values))
        print(f'Final equation for total value: {equation_str} = {total_value}')
        
    else:
        print('No solution found for the problem.')

# Run the solver
solve_multiple_knapsack()
<<<681>>>