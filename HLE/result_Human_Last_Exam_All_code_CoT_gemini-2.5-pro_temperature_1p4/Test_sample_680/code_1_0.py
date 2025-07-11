import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.
    """
    # --- Problem Definition ---
    city_width = 12  # km
    city_length = 11 # km
    total_budget = 45000  # usd

    # Tower B1 specs
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * (radius_b1 ** 2)

    # Tower B2 specs
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * (radius_b2 ** 2)

    # Total city area
    city_area = city_width * city_length

    # --- Brute-force Search for Optimal Solution ---
    best_n1 = 0
    best_n2 = 0
    max_coverage_area = 0

    # Determine maximum possible number of towers based on budget alone
    max_b1_possible = total_budget // cost_b1
    max_b2_possible = total_budget // cost_b2

    # Iterate through all feasible combinations of n1 and n2
    for n2 in range(max_b2_possible + 1):
        for n1 in range(max_b1_possible + 1):
            
            # Check budget constraint
            current_cost = (n1 * cost_b1) + (n2 * cost_b2)
            if current_cost > total_budget:
                continue

            # Calculate coverage area
            current_coverage_area = (n1 * area_b1) + (n2 * area_b2)

            # Check physical area constraint (non-overlapping)
            if current_coverage_area > city_area:
                continue
            
            # If this combination gives better coverage, update the best solution
            if current_coverage_area > max_coverage_area:
                max_coverage_area = current_coverage_area
                best_n1 = n1
                best_n2 = n2

    # --- Calculate and Print Results ---
    
    # Coverage ratio rounded to the nearest percentage
    coverage_ratio = (max_coverage_area / city_area) * 100
    coverage_percent = round(coverage_ratio)
    
    final_cost = best_n1 * cost_b1 + best_n2 * cost_b2

    print("Optimal Solution Found:")
    print(f"Number of B1 towers (b1): {best_n1}")
    print(f"Number of B2 towers (b2): {best_n2}\n")
    
    print("Final Equation Components:")
    # Using format() for cleaner alignment
    print("Cost Check:    {n1} * {c1} + {n2} * {c2} = {cost} <= {budget}".format(
        n1=best_n1, c1=cost_b1, n2=best_n2, c2=cost_b2, cost=final_cost, budget=total_budget))
    
    print("Coverage Area: {n1} * {a1:.2f} + {n2} * {a2:.2f} = {area:.2f} km^2".format(
        n1=best_n1, a1=area_b1, n2=best_n2, a2=area_b2, area=max_coverage_area))
        
    print("Coverage Ratio: ({area:.2f} km^2 / {city_a} km^2) * 100% = {c_ratio:.2f}%".format(
        area=max_coverage_area, city_a=city_area, c_ratio=coverage_ratio))
        
    print(f"\nRounded coverage (c): {coverage_percent}%\n")

    # Final answer in the required format
    print("Final Answer String:")
    print(f"{best_n1};{best_n2};{coverage_percent}")


solve_tower_placement()