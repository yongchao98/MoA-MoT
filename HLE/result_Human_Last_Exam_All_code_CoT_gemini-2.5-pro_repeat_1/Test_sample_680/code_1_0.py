import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # 1. Define problem parameters
    area_width = 12  # km
    area_length = 11 # km
    total_area = area_width * area_length

    # Tower B1 parameters
    cost1 = 1500
    radius1 = 1
    area1 = math.pi * radius1**2
    
    # Tower B2 parameters
    cost2 = 5000
    radius2 = 2
    area2 = math.pi * radius2**2

    # Budget constraint
    budget = 45000

    # 2. Solve the Integer Linear Program
    # We want to maximize the objective function: b1 * area1 + b2 * area2
    # This is equivalent to maximizing: b1 * 1 + b2 * 4
    # Subject to the constraint: b1 * 1500 + b2 * 5000 <= 45000
    # which simplifies to: 3 * b1 + 10 * b2 <= 90

    best_b1 = 0
    best_b2 = 0
    max_coverage_objective = -1

    # Iterate through all possible numbers of B2 towers
    max_b2 = int(budget / cost2)
    for b2 in range(max_b2 + 1):
        # For each number of B2 towers, calculate the maximum possible B1 towers
        remaining_budget = budget - (b2 * cost2)
        b1 = int(remaining_budget / cost1)
        
        # Calculate the value of the objective function (proportional to coverage)
        current_objective = b1 * (radius1**2) + b2 * (radius2**2)

        if current_objective > max_coverage_objective:
            max_coverage_objective = current_objective
            best_b1 = b1
            best_b2 = b2
    
    # 3. Calculate final coverage ratio
    max_coverage_area = best_b1 * area1 + best_b2 * area2
    coverage_ratio = max_coverage_area / total_area
    coverage_percentage = round(coverage_ratio * 100)

    # 4. Print the results as requested
    print("Optimization Problem Formulation:")
    print(f"Maximize Coverage = b1 * (pi * {radius1}^2) + b2 * (pi * {radius2}^2)")
    print(f"Subject to: b1 * {cost1} + b2 * {cost2} <= {budget}")
    print("\nOptimal Solution:")
    print(f"Number of B1 towers (b1) = {best_b1}")
    print(f"Number of B2 towers (b2) = {best_b2}")
    
    print("\nCoverage Ratio Calculation:")
    # Show the final equation with all numbers
    print(f"({best_b1} * pi * {radius1}^2 + {best_b2} * pi * {radius2}^2) / ({area_width} * {area_length}) = {coverage_ratio:.4f}")
    print(f"Coverage as a percentage = {coverage_percentage}%")
    
    print("\nFinal Answer:")
    # Print the answer in the required format
    print(f"{best_b1};{best_b2};{int(coverage_percentage)}")

solve_tower_placement()