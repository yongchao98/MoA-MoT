import math

def solve_chair_puzzle():
    """
    Calculates the maximum number of chairs that can be occupied
    based on the puzzle's rules.
    """
    # The total number of chairs in the row.
    n = 20

    # As explained in the plan, the maximum number of occupied chairs (k)
    # in a row of n chairs, with no two being adjacent, is given by the
    # formula k = floor((n + 1) / 2).

    # We apply this formula to find the solution.
    
    print(f"The problem is to find the maximum number of occupied chairs (k) out of a total of n = {n} chairs.")
    print("Based on the rules, the optimal arrangement has no two occupied chairs adjacent to each other.")
    print("\nThis can be solved using the formula: k <= (n + 1) / 2")
    
    print("\nLet's plug in the numbers:")
    # Show the substitution
    print(f"k <= ({n} + 1) / 2")
    
    # Show the calculation of the numerator
    numerator = n + 1
    print(f"k <= {numerator} / 2")

    # Show the result of the division
    division_result = numerator / 2
    print(f"k <= {division_result}")

    # Since the number of chairs must be a whole number, we take the floor.
    k = math.floor(division_result)
    
    print(f"As k must be an integer, the maximum value for k is {k}.")
    print("\n------------------------------------------------------")
    print(f"Final Answer: The maximum number of chairs that can be occupied is {k}.")
    print("------------------------------------------------------")


solve_chair_puzzle()