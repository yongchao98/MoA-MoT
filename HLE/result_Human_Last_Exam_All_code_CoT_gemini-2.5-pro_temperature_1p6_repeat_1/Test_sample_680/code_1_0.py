import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.

    The problem is formulated as an Integer Linear Program:
    Maximize Z = n1 + 4*n2
    Subject to:
    1. 3*n1 + 10*n2 <= 90 (Budget constraint)
    2. n1 + 4*n2 <= 33   (Physical packing constraint)
    where n1, n2 are non-negative integers.
    """
    
    # Total area of the city
    total_city_area = 12 * 11

    # Search range for n1 and n2 based on constraints
    # From budget: 3*n1 <= 90 => n1 <= 30. 10*n2 <= 90 => n2 <= 9
    max_n1 = 30
    max_n2 = 9
    
    best_n1 = 0
    best_n2 = 0
    max_objective_value = -1

    # Iterate through all possible values of n2 and n1
    for n2 in range(max_n2 + 1):
        for n1 in range(max_n1 + 1):
            # Check if the current combination (n1, n2) is feasible
            budget_ok = (3 * n1 + 10 * n2) <= 90
            packing_ok = (n1 + 4 * n2) <= 33
            
            if budget_ok and packing_ok:
                # Calculate the objective function value
                current_objective_value = n1 + 4 * n2
                
                # Check if this solution is better than the best one found so far.
                # Tie-breaking rule: if objective values are equal, prefer the one with more B2 towers.
                if current_objective_value > max_objective_value:
                    max_objective_value = current_objective_value
                    best_n1 = n1
                    best_n2 = n2
                elif current_objective_value == max_objective_value:
                    if n2 > best_n2:
                        best_n1 = n1
                        best_n2 = n2

    # Calculate the final coverage ratio
    b1 = best_n1
    b2 = best_n2
    
    # Total covered area = pi * (n1 * r1^2 + n2 * r2^2) = pi * (n1*1 + n2*4)
    total_coverage_area = math.pi * (b1 + 4 * b2)
    
    coverage_ratio = (total_coverage_area / total_city_area) * 100
    
    # c is the coverage ratio rounded to the nearest percentage
    c = int(round(coverage_ratio))
    
    # Print the result in the required format
    print(f"Optimal number of B1 towers (b1): {b1}")
    print(f"Optimal number of B2 towers (b2): {b2}")
    print(f"Maximized objective value (b1 + 4*b2): {b1 + 4 * b2}")
    print(f"Total cost: ${b1*1500 + b2*5000}")
    print(f"Total coverage area: {total_coverage_area:.2f} km^2")
    print(f"Total city area: {total_city_area} km^2")
    print(f"Coverage ratio: {coverage_ratio:.2f}%")
    print(f"\nFinal result in b1;b2;c format:")
    print(f"{b1};{b2};{c}")

solve_tower_placement()
<<<1;8;79>>>