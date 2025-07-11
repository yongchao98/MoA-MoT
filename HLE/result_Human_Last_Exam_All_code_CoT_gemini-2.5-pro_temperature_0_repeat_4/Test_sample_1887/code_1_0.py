def solve_set_theory_problem():
    """
    This function explains the solution to the set theory problem
    and prints the final order type.
    """

    print("Based on the problem's constraints, the set X contains all regular cardinals")
    print("Aleph_beta for which 1 <= beta < omega + 5.")
    print("")
    print("A cardinal Aleph_beta is regular if beta is a successor ordinal (or 0, or inaccessible).")
    print("In the range 1 <= beta < omega + 5, this means beta must be a successor ordinal.")
    print("")
    print("The successor ordinals beta < omega + 5 are:")
    print("1. The finite successors: 1, 2, 3, ...")
    print("2. The transfinite successors: omega+1, omega+2, omega+3, omega+4")
    print("")
    print("This divides the set X into two ordered parts:")
    part1_description = "The set {Aleph_1, Aleph_2, Aleph_3, ...}, which has order type omega."
    part2_description = "The set {Aleph_{omega+1}, Aleph_{omega+2}, Aleph_{omega+3}, Aleph_{omega+4}}, which has order type 4."
    
    print("Part 1:", part1_description)
    print("Part 2:", part2_description)
    print("")
    
    print("The order type of X is the ordinal sum of the order types of these two parts.")
    
    # The components of the final equation
    part1_type = "omega"
    operator = "+"
    part2_type = "4"
    
    print("\nFinal Equation:")
    print(f"{part1_type} {operator} {part2_type}")

solve_set_theory_problem()