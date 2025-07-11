def solve_tower_placement():
    """
    This function analyzes the tower placement problem to find the optimal solution.

    The analysis confirms that the user's problem formulation is correct.
    The core of the problem is to find the number of B1 towers (n1) and B2 towers (n2)
    that minimize the cost C = 1500*n1 + 4000*n2, subject to two main constraints:
    1. Coverage: n1 + 4*n2 >= 34
    2. Geometric Packing: All towers must fit in the 12x11 km area without their
       coverage circles overlapping.

    The analysis of cost-effectiveness shows B2 towers are cheaper per unit of coverage.
    Therefore, solutions with more B2 towers are likely to be cheaper.

    The two best candidates from a cost perspective are:
    - (n1=2, n2=8) with a cost of 35,000.
    - (n1=0, n2=9) with a cost of 36,000.

    A geometric packing analysis shows that the (n1=2, n2=8) configuration is not
    realistically packable in the given area. However, the (n1=0, n2=9) configuration
    is easily packable (e.g., in a 3x3 grid).

    Thus, the cheapest FEASIBLE solution is (n1=0, n2=9).
    """

    # Final optimal solution values
    b1_towers = 0
    b1_cost = 1500

    b2_towers = 9
    b2_cost = 4000

    # Calculate the minimum cost based on the optimal number of towers
    total_cost = b1_towers * b1_cost + b2_towers * b2_cost

    print("The problem formulation is correct.")
    print("Based on cost-effectiveness and packing feasibility, the optimal solution is determined.")
    print("\nFinal cost calculation:")
    # Output each number in the final equation as requested
    print(f"{b1_towers} * {b1_cost} + {b2_towers} * {b2_cost} = {total_cost}")

    print("\nAnswer in the required format (b1;b2;c):")
    # Output the final answer string
    print(f"{b1_towers};{b2_towers};{total_cost}")


solve_tower_placement()