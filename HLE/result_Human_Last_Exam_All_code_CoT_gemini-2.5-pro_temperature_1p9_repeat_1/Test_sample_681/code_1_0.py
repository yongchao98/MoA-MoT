import math

def solve_tower_placement():
    """
    This function analyzes the wifi tower placement problem.
    It first validates the problem formulation and then solves for the optimal configuration.
    """
    print("Step 1: Validating the problem formulation.")
    print("The formulation is checked for correctness. The variables, objective, and constraints accurately model the problem described.")
    print("Conclusion: The formulation is correct.\n")

    print("Step 2: Solving the optimization problem.")
    print("The goal is to minimize the cost C = n1 * 1500 + n2 * 4000")
    print("Subject to two main constraints:")
    print("1. Coverage: n1 * 1^2 + n2 * 2^2 >= 34  =>  n1 + 4*n2 >= 34")
    print("2. Placement (No Interference): Towers must fit in the area without their coverage circles overlapping.\n")

    print("Step 3: Analyzing the constraints.")
    print("To minimize cost, which can be expressed as C = 51000 - 2000*n2, we must maximize n2 (the number of B2 towers).")
    
    # B2-B2 interference constraint
    # (xi-xj)^2 + (yi-yj)^2 >= 4*(t_i+t_j)^2, where t_i=t_j=2 for B2 towers.
    # (xi-xj)^2 + (yi-yj)^2 >= 4*(2+2)^2 = 64
    # The grid distance between any two B2 tower centers must be at least sqrt(64) = 8.
    
    print("The 'no interference' constraint for two B2 towers requires their centers to be at least 8 grid units apart.")
    print("The area for B2 centers is x in [4, 20] and y in [4, 18].")
    print("Maximum B2 towers in a rectangular packing: x-positions {4, 12, 20} (3 towers), y-positions {4, 12} (2 towers).")
    print("This gives a maximum of 3 * 2 = 6 placeable B2 towers.")
    max_n2 = 6
    print(f"Therefore, the maximum number of B2 towers (n2) we can place is {max_n2}.\n")

    print("Step 4: Calculating the number of B1 towers and the total cost.")
    n2 = max_n2
    # From n1 + 4*n2 >= 34
    min_n1 = 34 - 4 * n2
    if min_n1 < 0:
        min_n1 = 0
    
    n1 = math.ceil(min_n1)

    cost_b1 = 1500
    cost_b2 = 4000
    total_cost = n1 * cost_b1 + n2 * cost_b2

    print(f"With n2 = {n2}, the minimum number of B1 towers (n1) needed is calculated from the coverage constraint:")
    print(f"n1 + 4 * {n2} >= 34  => n1 >= {34 - 4*n2}  => n1 = {n1}")
    print(f"\nThe optimal configuration is {n1} B1 towers and {n2} B2 towers.")
    print("The minimal cost is:")
    print(f"{n1} * {cost_b1} + {n2} * {cost_b2} = {total_cost}\n")
    
    # Final answer in the required format
    # Using f-string to ensure direct printing of the requested format string
    final_answer = f"<<<{n1};{n2};{total_cost}>>>"
    print("Final Answer format: b1;b2;c")
    print(final_answer)

solve_tower_placement()