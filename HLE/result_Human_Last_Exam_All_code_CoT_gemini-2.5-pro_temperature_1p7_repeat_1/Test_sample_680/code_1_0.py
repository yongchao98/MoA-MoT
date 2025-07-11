import math

def solve_tower_optimization():
    """
    Solves the WiFi tower placement optimization problem.
    """
    # --- Problem Parameters ---
    city_area = 12 * 11  # km^2
    budget_total = 45000

    # B1 Tower parameters
    radius1 = 1.0
    cost1 = 1500.0
    area1 = math.pi * radius1**2

    # B2 Tower parameters
    radius2 = 2.0
    cost2 = 5000.0
    area2 = math.pi * radius2**2

    # --- ILP Formulation & Solver ---
    # Variables: n1, n2 (number of towers)
    # Maximize: Z = n1 + 4*n2
    # Constraints:
    # 1. 3*n1 + 10*n2 <= 90 (Budget)
    # 2. n1 + 4*n2 <= 42 (Area)
    
    best_n1 = 0
    best_n2 = 0
    max_objective_value = -1

    # Iterate through all possible numbers of B2 towers (n2)
    # From budget constraint: 10*n2 <= 90 => n2 <= 9
    # From area constraint: 4*n2 <= 42 => n2 <= 10
    # Thus, the maximum possible value for n2 is 9.
    for n2 in range(10): # from 0 to 9
        
        # Calculate max n1 allowed by budget constraint
        # 3*n1 <= 90 - 10*n2
        max_n1_budget = (90 - 10 * n2) / 3
        
        # Calculate max n1 allowed by area constraint
        # n1 <= 42 - 4*n2
        max_n1_area = 42 - 4 * n2

        # n1 must be a non-negative integer satisfying both constraints
        if max_n1_budget < 0 or max_n1_area < 0:
            continue
            
        n1 = math.floor(min(max_n1_budget, max_n1_area))

        # Calculate the objective function value
        current_objective_value = n1 + 4 * n2
        
        # If this solution is better than the previous best, update it
        if current_objective_value > max_objective_value:
            max_objective_value = current_objective_value
            best_n1 = n1
            best_n2 = n2

    # --- Calculate and Print Final Results ---
    final_cost = cost1 * best_n1 + cost2 * best_n2
    final_coverage_area = area1 * best_n1 + area2 * best_n2
    coverage_ratio = (final_coverage_area / city_area) * 100
    
    # Rounded to the nearest percentage
    coverage_percentage = round(coverage_ratio)

    print("--- Optimal Solution Found ---")
    print(f"Number of B1 towers (b1): {best_n1}")
    print(f"Number of B2 towers (b2): {best_n2}\n")

    print("--- Verifying Constraints ---")
    # Budget Equation: 1500 * n1 + 5000 * n2 <= 45000
    print(f"Total Cost: {cost1:.0f}*{best_n1} + {cost2:.0f}*{best_n2} = {final_cost:.0f} (<= {budget_total})")
    # Area Equation: π*(1^2*n1 + 2^2*n2) <= 132
    print(f"Total Coverage: π*({radius1**2:.0f}²*{best_n1} + {radius2**2:.0f}²*{best_n2}) = {final_coverage_area:.2f} m² (<= {city_area})\n")

    print("--- Final Coverage Ratio ---")
    # Coverage Ratio Equation: (Total Coverage / City Area) * 100
    print(f"Ratio: ({final_coverage_area:.2f} / {city_area}) * 100 = {coverage_ratio:.2f}%")
    print(f"Rounded to nearest percentage (c): {coverage_percentage}%\n")

    # Output the final answer in the required format: b1;b2;c
    final_answer = f"{best_n1};{best_n2};{coverage_percentage}"
    print(f"Final Formatted Answer: {final_answer}")
    
    return final_answer

final_answer_string = solve_tower_optimization()
print(f'<<<{final_answer_string}>>>')
