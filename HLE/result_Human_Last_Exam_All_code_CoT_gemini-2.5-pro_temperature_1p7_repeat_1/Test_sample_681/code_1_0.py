import math

def solve():
    """
    Finds the optimal number of B1 and B2 towers to minimize cost
    subject to a total coverage constraint.
    """
    cost_b1 = 1500
    cost_b2 = 4000
    coverage_sq_b1 = 1**2
    coverage_sq_b2 = 2**2
    
    total_area = 12 * 11
    required_coverage_ratio = 0.80
    required_total_coverage = total_area * required_coverage_ratio
    # In the problem formulation, this is Σ(t_i^2) >= 105.6 / pi
    required_coverage_units_sq = math.ceil(required_total_coverage / math.pi)

    min_cost = float('inf')
    optimal_n1 = -1
    optimal_n2 = -1

    # Iterate through a reasonable range of B2 towers
    # Max n2 is ceil(required_coverage_units_sq / coverage_sq_b2) + a small margin
    max_n2_check = math.ceil(required_coverage_units_sq / coverage_sq_b2) + 2
    
    for n2 in range(int(max_n2_check)):
        # Calculate coverage provided by n2 B2 towers
        coverage_from_b2 = n2 * coverage_sq_b2
        
        # Calculate remaining coverage needed from B1 towers
        remaining_coverage_needed = required_coverage_units_sq - coverage_from_b2
        
        n1 = 0
        if remaining_coverage_needed > 0:
            # Each B1 tower provides coverage_sq_b1 (1 unit)
            n1 = math.ceil(remaining_coverage_needed / coverage_sq_b1)

        # Calculate the cost for this combination
        current_cost = n1 * cost_b1 + n2 * cost_b2

        # If this cost is the new minimum, store this combination
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_n1 = n1
            optimal_n2 = n2

    print("Step 1: Determine the most cost-effective combination of towers based on coverage.")
    print(f"Coverage requirement Σ(t_i^2) >= {required_coverage_units_sq}")
    print(f"Optimal combination found: {optimal_n1} B1 towers and {optimal_n2} B2 towers.")
    
    print("\nStep 2: Verify the constraints for this combination.")
    print("Coverage constraint check:")
    total_coverage_units = optimal_n1 * coverage_sq_b1 + optimal_n2 * coverage_sq_b2
    print(f"Equation: n1*{coverage_sq_b1} + n2*{coverage_sq_b2} >= {required_coverage_units_sq}")
    print(f"Values:   {optimal_n1}*{coverage_sq_b1} + {optimal_n2}*{coverage_sq_b2} = {total_coverage_units} (which is >= {required_coverage_units_sq})")
    
    print("\nObjective function (cost) calculation:")
    print(f"Equation: Cost = n1*{cost_b1} + n2*{cost_b2}")
    print(f"Values:   Cost = {optimal_n1}*{cost_b1} + {optimal_n2}*{cost_b2} = {min_cost}")

    print("\nStep 3: Confirm placement feasibility (as discussed in the text, this combination is packable).")

    # Final result in the required format.
    print("\n---")
    print("Final Answer in requested format b1;b2;c:")
    print(f"{optimal_n1};{optimal_n2};{min_cost}")

solve()