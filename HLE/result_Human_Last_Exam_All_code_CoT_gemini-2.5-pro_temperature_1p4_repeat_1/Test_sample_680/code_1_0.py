import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage area.
    """
    # --- Problem Definition ---
    # City area
    rect_w = 12  # km
    rect_h = 11  # km
    rect_area = rect_w * rect_h

    # Total budget
    total_budget = 45000  # usd

    # B1 Tower specifications
    b1_radius = 1
    b1_cost = 1500
    b1_area = math.pi * b1_radius**2

    # B2 Tower specifications
    b2_radius = 2
    b2_cost = 5000
    b2_area = math.pi * b2_radius**2

    # --- Optimization ---
    # We will use a greedy approach based on the area-per-cost efficiency.
    # B2 tower efficiency: 4*pi / 5000 = 0.0008*pi area/usd
    # B1 tower efficiency: pi / 1500 approx 0.00067*pi area/usd
    # B2 is more efficient, so we prioritize it.

    best_b1 = 0
    best_b2 = 0
    max_coverage_area = 0
    
    # Iterate b2 from the max possible downwards
    max_b2_possible = total_budget // b2_cost
    for num_b2 in range(max_b2_possible, -1, -1):
        # Calculate budget remaining after purchasing b2 towers
        remaining_budget = total_budget - (num_b2 * b2_cost)
        
        # Calculate how many b1 towers can be purchased with the rest
        num_b1 = remaining_budget // b1_cost
        
        # Calculate the coverage area for this combination
        current_coverage_area = (num_b1 * b1_area) + (num_b2 * b2_area)
        
        # Check simplified non-overlap (area) constraint
        if current_coverage_area <= rect_area:
            # Since we iterate from the most efficient tower option downwards,
            # the first combination that satisfies the constraints is the optimal one.
            best_b1 = num_b1
            best_b2 = num_b2
            max_coverage_area = current_coverage_area
            break # Exit loop once the best solution is found
            
    # --- Calculate Final Results ---
    final_cost = best_b1 * b1_cost + best_b2 * b2_cost
    coverage_ratio = max_coverage_area / rect_area
    coverage_percentage = round(coverage_ratio * 100)

    # --- Output Results ---
    print("Optimization Problem Formulation:")
    print(f"Maximize: pi * (b1 + 4*b2)")
    print(f"Subject to:")
    print(f"  {b1_cost}*b1 + {b2_cost}*b2 <= {total_budget}  (Budget)")
    print(f"  pi*(b1 + 4*b2) <= {rect_area}            (Area)\n")

    print("Optimal Solution Found:")
    print("Final Equation for Cost:")
    print(f"{best_b1} * ${b1_cost} + {best_b2} * ${b2_cost} = ${final_cost}")
    
    print("\nFinal Equation for Coverage:")
    print(f"({best_b1} * {b1_area:.2f}km^2 + {best_b2} * {b2_area:.2f}km^2) / {rect_area}km^2 = {coverage_ratio:.2%}")

    print("\nFinal Answer:")
    print(f"{best_b1};{best_b2};{coverage_percentage}")


solve_tower_placement()