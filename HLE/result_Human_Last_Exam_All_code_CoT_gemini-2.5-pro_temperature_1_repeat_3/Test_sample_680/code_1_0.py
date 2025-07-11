import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # Problem parameters
    rect_area = 12 * 11  # km^2
    budget = 45000       # usd

    # Tower B1 parameters
    radius1 = 1          # km
    cost1 = 1500         # usd
    area1 = math.pi * (radius1 ** 2)

    # Tower B2 parameters
    radius2 = 2          # km
    cost2 = 5000         # usd
    area2 = math.pi * (radius2 ** 2)

    # Variables to store the best solution found
    best_b1 = 0
    best_b2 = 0
    max_coverage = 0.0

    # Determine the maximum possible number of each tower type to set loop bounds
    max_b1_possible = budget // cost1
    max_b2_possible = budget // cost2

    # Iterate through all possible numbers of towers
    # We prioritize b2 as it provides more coverage per tower
    for b2 in range(max_b2_possible + 1):
        for b1 in range(max_b1_possible + 1):
            
            # Check budget constraint
            current_cost = b1 * cost1 + b2 * cost2
            if current_cost > budget:
                # Since the inner loop for b1 is increasing, we can break
                # early if the budget is exceeded.
                break

            # Check non-overlapping area constraint
            total_circle_area = b1 * area1 + b2 * area2
            if total_circle_area > rect_area:
                continue

            # Check if this combination gives better coverage
            if total_circle_area > max_coverage:
                max_coverage = total_circle_area
                best_b1 = b1
                best_b2 = b2
    
    # Calculate the final coverage ratio
    coverage_ratio = (max_coverage / rect_area) * 100
    
    # Output the results and the final equations
    print("Optimization Problem Formulation:")
    print("Maximize: b1 * pi * (1^2) + b2 * pi * (2^2)")
    print("Subject to:")
    print(f"  1500*b1 + 5000*b2 <= {budget}  (Budget Constraint)")
    print(f"  b1*pi*(1^2) + b2*pi*(2^2) <= {rect_area}  (Area Constraint)")
    print("\n--- Optimal Solution Found ---\n")
    print(f"Number of B1 towers (b1): {best_b1}")
    print(f"Number of B2 towers (b2): {best_b2}")
    
    print("\n--- Final Equation Values ---")
    final_cost = best_b1 * cost1 + best_b2 * cost2
    print(f"Total Cost = {best_b1} * {cost1} + {best_b2} * {cost2} = ${final_cost}")
    print(f"Total Coverage = {best_b1} * pi*{radius1}^2 + {best_b2} * pi*{radius2}^2 = {max_coverage:.2f} km^2")
    
    coverage_percentage = int(round(coverage_ratio))
    print(f"\nCoverage Ratio (c) = ({max_coverage:.2f} / {rect_area}) * 100% = {coverage_ratio:.2f}%")
    print(f"Rounded to nearest percentage: {coverage_percentage}%\n")

    # Final answer in the required format
    print("Final Answer Format (b1;b2;c):")
    print(f"{best_b1};{best_b2};{coverage_percentage}")


solve_tower_placement()
<<<0;9;86>>>