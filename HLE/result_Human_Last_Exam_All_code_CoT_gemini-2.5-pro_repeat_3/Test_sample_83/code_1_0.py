from ortools.linear_solver import pywraplp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using Google's OR-Tools library.
    """
    # 1. Define the data
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(capacities)

    # 2. Create the solver
    # Using the SCIP backend, which is good for mixed-integer programming
    solver = pywraplp.Solver.CreateSolver('SCIP')
    if not solver:
        print("SCIP solver not available.")
        return

    # 3. Create the variables
    # x[i, j] is a binary variable, 1 if item i is put in knapsack j.
    x = {}
    for i in range(num_items):
        for j in range(num_knapsacks):
            x[i, j] = solver.IntVar(0, 1, f'x_{i}_{j}')

    # 4. Define the constraints
    # a) Each item can be assigned to at most one knapsack.
    for i in range(num_items):
        solver.Add(sum(x[i, j] for j in range(num_knapsacks)) <= 1)

    # b) The total weight in each knapsack cannot exceed its capacity.
    for j in range(num_knapsacks):
        solver.Add(
            sum(weights[i] * x[i, j] for i in range(num_items)) <= capacities[j]
        )

    # 5. Define the objective function
    # We want to maximize the total value of the items in the knapsacks.
    objective = solver.Objective()
    for i in range(num_items):
        for j in range(num_knapsacks):
            objective.SetCoefficient(x[i, j], values[i])
    objective.SetMaximization()

    # 6. Solve the problem
    status = solver.Solve()

    # 7. Print the solution
    if status == pywraplp.Solver.OPTIMAL:
        total_value = int(solver.Objective().Value())
        print(f"Maximum possible total value: {total_value}")

        selected_values = []
        total_weight_knapsacks = [0] * num_knapsacks

        print("\n--- Solution Details ---")
        for j in range(num_knapsacks):
            knapsack_items = []
            knapsack_value = 0
            knapsack_weight = 0
            for i in range(num_items):
                if x[i, j].solution_value() > 0.5: # Item i is in knapsack j
                    selected_values.append(values[i])
                    knapsack_items.append(f"Item {i+1} (Value: {values[i]}, Weight: {weights[i]})")
                    knapsack_value += values[i]
                    knapsack_weight += weights[i]

            print(f"\nKnapsack {j+1} (Capacity: {capacities[j]}):")
            if knapsack_items:
                for item_desc in knapsack_items:
                    print(f"  - {item_desc}")
            else:
                print("  - Empty")
            print(f"  Total Value: {knapsack_value}, Total Weight: {knapsack_weight}")
        
        # Sort the selected values for a clean equation
        selected_values.sort(reverse=True)
        equation_str = " + ".join(map(str, selected_values))

        print("\n--- Final Equation ---")
        print(f"{equation_str} = {total_value}")

    else:
        print("The problem does not have an optimal solution.")

if __name__ == '__main__':
    solve_multiple_knapsack()
<<<694>>>