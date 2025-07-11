def solve_arboricity_bounds():
    """
    This function determines the complexity classes for the arboricity bounds f1(n) and f2(n).

    f1(n) corresponds to the case c=1. The analysis shows that the arboricity
    is bounded by Theta(log(n)/log(log(n))). This complexity class is
    omega(sqrt(log(n))) but o(log(n)), which corresponds to category 4.

    f2(n) corresponds to the case c=2. The analysis shows that the arboricity
    is bounded by O(1) with high probability. This corresponds to category 1.
    
    The resulting two-digit number is the concatenation of the categories for f1 and f2.
    """

    # Category for f1(n) where c = 1
    f1_category = 4

    # Category for f2(n) where c = 2
    f2_category = 1

    # The problem asks for a two-digit number.
    # The first digit is for f1, the second is for f2.
    final_answer = int(f"{f1_category}{f2_category}")
    
    print(final_answer)

solve_arboricity_bounds()