import math

def solve_tower_placement():
    """
    Solves the Wi-Fi tower placement optimization problem.
    """
    # --- Problem Definition ---
    # Area dimensions
    rect_width = 12
    rect_length = 11
    total_area = rect_width * rect_length

    # Tower B1 properties
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * (radius_b1 ** 2)
    # Bounding box area for placement constraint
    bbox_area_b1 = (2 * radius_b1) ** 2

    # Tower B2 properties
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * (radius_b2 ** 2)
    # Bounding box area for placement constraint
    bbox_area_b2 = (2 * radius_b2) ** 2

    # Budget constraint
    total_budget = 45000

    # --- Optimization ---
    best_b1 = 0
    best_b2 = 0
    max_coverage_score = -1
    min_cost_at_max_coverage = float('inf')

    # Determine the search range for b1 and b2
    max_b1 = min(total_budget // cost_b1, total_area // bbox_area_b1)
    max_b2 = min(total_budget // cost_b2, total_area // bbox_area_b2)

    # Iterate through all possible combinations of towers
    for b1 in range(max_b1 + 1):
        for b2 in range(max_b2 + 1):
            
            # Calculate total cost and total bounding box area
            current_cost = b1 * cost_b1 + b2 * cost_b2
            current_bbox_area = b1 * bbox_area_b1 + b2 * bbox_area_b2

            # Check if the combination is valid (within budget and area constraints)
            if current_cost <= total_budget and current_bbox_area <= total_area:
                
                # The objective is to maximize (b1 * pi + b2 * 4*pi), which is equivalent
                # to maximizing the coverage score (b1 + 4*b2).
                coverage_score = b1 + 4 * b2
                
                # Check if this solution is better than the current best
                if coverage_score > max_coverage_score:
                    max_coverage_score = coverage_score
                    min_cost_at_max_coverage = current_cost
                    best_b1 = b1
                    best_b2 = b2
                # If coverage is the same, check if it's cheaper (tie-breaker)
                elif coverage_score == max_coverage_score:
                    if current_cost < min_cost_at_max_coverage:
                        min_cost_at_max_coverage = current_cost
                        best_b1 = b1
                        best_b2 = b2
    
    # --- Final Calculation ---
    # Total coverage area from the optimal combination
    total_coverage_area = best_b1 * area_b1 + best_b2 * area_b2
    
    # Coverage ratio as a percentage
    coverage_ratio = (total_coverage_area / total_area) * 100
    
    # Round to the nearest percentage
    coverage_percentage = round(coverage_ratio)

    # Print the result in the specified format: b1;b2;c
    print(f"{best_b1};{best_b2};{coverage_percentage}")

solve_tower_placement()