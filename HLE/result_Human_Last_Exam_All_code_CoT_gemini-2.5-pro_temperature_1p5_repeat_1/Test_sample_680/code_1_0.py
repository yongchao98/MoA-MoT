import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.
    """
    # City and Budget parameters
    city_area = 12 * 11
    total_budget = 45000

    # Tower parameters
    b1_radius = 1
    b1_cost = 1500
    b1_area = math.pi * b1_radius**2
    # The "value" is proportional to area, normalized by pi
    b1_value_factor = b1_radius**2

    b2_radius = 2
    b2_cost = 5000
    b2_area = math.pi * b2_radius**2
    b2_value_factor = b2_radius**2

    # Packing density heuristic (e.g., 80%)
    max_packing_density = 0.8
    max_total_circle_area = city_area * max_packing_density

    best_b1 = 0
    best_b2 = 0
    max_value = -1

    # Iterate through all possible numbers of the more expensive tower (B2)
    max_b2_by_budget = total_budget // b2_cost
    for n2 in range(max_b2_by_budget + 1):
        # Calculate remaining budget after placing n2 B2 towers
        remaining_budget = total_budget - n2 * b2_cost
        
        # Calculate max B1 towers possible with remaining budget
        n1_by_budget = remaining_budget // b1_cost
        
        # Calculate max B1 towers possible based on packing density
        remaining_area_for_packing = max_total_circle_area - n2 * b2_area
        if remaining_area_for_packing < 0:
            n1_by_packing = -1 # Impossible to pack
        else:
            n1_by_packing = math.floor(remaining_area_for_packing / b1_area)
            
        # The number of B1 towers is limited by the stricter of the two constraints
        n1 = min(n1_by_budget, n1_by_packing)

        if n1 < 0:
            continue

        # Calculate the objective value to maximize
        current_value = n1 * b1_value_factor + n2 * b2_value_factor
        
        # If this combination is better, store it
        if current_value > max_value:
            max_value = current_value
            best_b1 = n1
            best_b2 = n2

    # Calculate the final coverage percentage
    total_coverage_area = best_b1 * b1_area + best_b2 * b2_area
    coverage_ratio = (total_coverage_area / city_area) * 100
    coverage_percentage = round(coverage_ratio)
    
    # Output the results
    print(f"Optimal number of B1 towers (b1): {best_b1}")
    print(f"Optimal number of B2 towers (b2): {best_b2}")
    print(f"Maximized coverage value: {best_b1} * {int(b1_value_factor)} + {best_b2} * {int(b2_value_factor)} = {int(max_value)}")
    print(f"Total coverage area: {total_coverage_area:.2f} km^2")
    print(f"City area: {city_area} km^2")
    print(f"Coverage Ratio: {coverage_ratio:.2f}%")
    print("---")
    print("Final Answer:")
    # Final answer in the format b1;b2;c
    print(f"{best_b1};{best_b2};{coverage_percentage}")
    
solve_tower_placement()