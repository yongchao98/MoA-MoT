import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # 1. Define Variables and Constants
    city_width = 12  # km
    city_height = 11  # km
    total_area = city_width * city_height
    budget = 45000  # usd

    # Tower B1 parameters
    cost1 = 1500
    radius1 = 1
    area1 = math.pi * radius1**2

    # Tower B2 parameters
    cost2 = 5000
    radius2 = 2
    area2 = math.pi * radius2**2

    # 2. Print the Formulation
    print("This script formulates and solves the tower placement optimization problem.")
    print("-" * 30)
    print("Problem Formulation:")
    print("Let b1 = number of B1 towers, and b2 = number of B2 towers.")
    print("\nObjective: Maximize Total Coverage")
    # Using f-string to display the equation with calculated values
    print(f"Maximize: b1 * {area1:.2f} + b2 * {area2:.2f}")

    print("\nSubject to Constraints:")
    print(f"1. Budget: b1 * {cost1} + b2 * {cost2} <= {budget}")
    print(f"2. Area:   b1 * {area1:.2f} + b2 * {area2:.2f} <= {total_area}")
    print("-" * 30)

    # 3. Solve the Integer Program by iterating through feasible solutions
    best_b1 = 0
    best_b2 = 0
    max_coverage = 0

    # Determine the maximum possible number of b2 towers by budget
    max_b2_by_budget = budget // cost2

    # Iterate through all possible numbers of b2 towers
    for b2 in range(max_b2_by_budget + 1):
        # For the current number of b2 towers, find the max b1 towers we can afford
        remaining_budget = budget - (b2 * cost2)
        b1 = remaining_budget // cost1

        # Calculate the total cost and coverage for this combination
        current_cost = b1 * cost1 + b2 * cost2
        current_coverage = b1 * area1 + b2 * area2

        # Check if the area constraint is met and if this is the best solution so far
        if current_coverage <= total_area:
            if current_coverage > max_coverage:
                max_coverage = current_coverage
                best_b1 = b1
                best_b2 = b2

    # 4. Calculate Final Ratio
    coverage_ratio = (max_coverage / total_area)
    coverage_percentage = round(coverage_ratio * 100)

    # 5. Display the Result
    print("Optimal Solution Found:")
    print(f"Number of B1 towers (b1) = {best_b1}")
    print(f"Number of B2 towers (b2) = {best_b2}")
    
    final_cost = best_b1 * cost1 + best_b2 * cost2
    print(f"Total Cost = {best_b1} * {cost1} + {best_b2} * {cost2} = ${final_cost}")
    print("-" * 30)
    
    print("Final Coverage Calculation:")
    print(f"Maximized Coverage = {best_b1} * {area1:.2f} + {best_b2} * {area2:.2f} = {max_coverage:.2f} km^2")
    print(f"Total City Area = {total_area} km^2")
    print(f"Coverage Ratio = {max_coverage:.2f} / {total_area} = {coverage_ratio:.4f}")
    print(f"Coverage Percentage (c) = {coverage_ratio:.4f} * 100 = {round(coverage_ratio * 100, 2)}%")
    print(f"Rounded to nearest percentage = {coverage_percentage}%")
    print("-" * 30)
    
    # Final answer in the required format
    final_answer = f"{best_b1};{best_b2};{coverage_percentage}"
    print(f"Final Answer (b1;b2;c): {final_answer}")


if __name__ == '__main__':
    solve_tower_optimization()