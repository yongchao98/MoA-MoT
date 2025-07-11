import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.
    """
    # Problem constants
    rect_width = 12
    rect_height = 11
    total_area = rect_width * rect_height
    budget = 45000

    # B1 Tower specifications
    cost1 = 1500
    radius1 = 1
    
    # B2 Tower specifications
    cost2 = 5000
    radius2 = 2

    # Search for the optimal number of towers
    best_n1 = 0
    best_n2 = 0
    max_weighted_coverage = -1
    min_cost_for_max_coverage = float('inf')

    # Iterate through possible numbers of B2 towers
    # The maximum possible n2 is limited by the budget (45000 / 5000 = 9)
    max_b2_towers = budget // cost2
    for n2 in range(max_b2_towers, -1, -1):
        
        # Given n2, find the maximum possible n1 that satisfies both constraints
        
        # 1. Budget constraint: 3*n1 + 10*n2 <= 90 => n1 <= (90 - 10*n2) / 3
        n1_from_budget = (90 - 10 * n2) // 3
        
        # 2. Packing constraint: n1 + 4*n2 <= 33 => n1 <= 33 - 4*n2
        n1_from_packing = 33 - 4 * n2

        # To be feasible, n1 must be positive and satisfy both constraints
        if n1_from_budget < 0 or n1_from_packing < 0:
            continue
            
        n1 = min(n1_from_budget, n1_from_packing)
        
        current_weighted_coverage = n1 + 4 * n2
        current_cost = n1 * cost1 + n2 * cost2

        # Check if this solution is better than the current best
        if current_weighted_coverage > max_weighted_coverage:
            max_weighted_coverage = current_weighted_coverage
            best_n1 = n1
            best_n2 = n2
            min_cost_for_max_coverage = current_cost
        # Tie-breaking rule: if coverage is the same, choose the cheaper option
        elif current_weighted_coverage == max_weighted_coverage:
            if current_cost < min_cost_for_max_coverage:
                best_n1 = n1
                best_n2 = n2
                min_cost_for_max_coverage = current_cost
                
    # Calculate the final coverage ratio
    total_coverage_area = (best_n1 * math.pi * radius1**2) + (best_n2 * math.pi * radius2**2)
    coverage_ratio = total_coverage_area / total_area
    coverage_percentage = round(coverage_ratio * 100)
    
    # Output the equation details as requested
    pi_str = f"{math.pi:.2f}"
    print(f"The integer linear programming formulation is:")
    print(f"Maximize: n1 * {pi_str} * {radius1}^2 + n2 * {pi_str} * {radius2}^2")
    print(f"Subject to:")
    print(f"  {cost1}*n1 + {cost2}*n2 <= {budget}  (Budget)")
    print(f"  n1 + 4*n2 <= 33              (Simplified Packing)")
    print("\nSolving this yields the optimal numbers of towers.")
    
    print("\nCalculation of the final coverage percentage:")
    print(f"Coverage Ratio = (({best_n1} * {pi_str} * {radius1}**2) + ({best_n2} * {pi_str} * {radius2}**2)) / ({rect_width} * {rect_height})")
    print(f"Result = {coverage_ratio:.4f} or {coverage_percentage}%\n")

    # Final answer in the specified format
    final_answer = f"{best_n1};{best_n2};{coverage_percentage}"
    print(f"Final Answer (b1;b2;c): {final_answer}")
    
    # Required final output format
    return final_answer

# Execute the function and print the final answer block
final_answer_str = solve_tower_placement()
print(f"<<<{final_answer_str}>>>")