import math

def solve():
    """
    This function solves the problem based on the derived formula.
    """
    
    # Part (a): Determine the expression for the maximum number of roots.
    # The mathematical derivation shows the maximum number of roots of R_t in ]0, 1[
    # is given by the formula t*(t-1)/2.
    expression_a = "t*(t-1)/2"

    # Part (b): Calculate the maximum number of roots for t = 5.
    t = 5
    
    # Applying the formula from part (a).
    # The problem asks to output the numbers in the final equation. We will show the calculation in comments.
    # Calculation: 5 * (5 - 1) / 2 = 5 * 4 / 2 = 20 / 2 = 10.
    result_b = t * (t - 1) // 2

    # Print the answer in the specified format.
    print(f"(a) {expression_a}; (b) {result_b}")

solve()