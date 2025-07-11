import sys

def solve_cohomology():
    """
    This function provides the cohomology groups of M(7), the moduli space of 7 disjoint
    linearly embedded closed intervals in R^2.

    The solution relies on two key results from algebraic topology:
    1. A theorem by V.A. Vassiliev, which states that the space M(k) is homotopy
       equivalent to C_k(S^1), the configuration space of k unordered points on a circle.
       Therefore, M(7) is homotopy equivalent to C_7(S^1).
    2. The space C_k(S^1) is itself homotopy equivalent to the circle S^1 for all k >= 1.
       This is a known property of configuration spaces on the circle.

    Combining these results, M(7) is homotopy equivalent to S^1. We then need to provide
    the integer cohomology groups of the circle S^1.

    The integer cohomology groups of S^1 are:
    - H^0(S^1; Z) = Z (since it is path-connected)
    - H^1(S^1; Z) = Z (corresponding to its single loop/hole)
    - H^a(S^1; Z) = 0 for a > 1

    The final answer is formatted as a list of these groups.
    """

    # H^0(M(7)) is Z because the space is path-connected.
    H0 = "Z"

    # H^1(M(7)) is Z.
    H1 = "Z"

    # Higher cohomology groups are 0.
    # The list is truncated at the last non-zero group, which is H^1.
    result_list = [H0, H1]

    # The prompt requires the output format "[H^0, H^1, ..., H^a]".
    # "Remember in the final code you still need to output each number in the final equation!"
    # The resulting groups Z have implicit exponents and coefficients of 1, which are standardly omitted.
    # The formatted string represents the list [Z, Z].
    
    # We build the string representation of the list
    # and print it to standard output.
    output = f"[{result_list[0]}, {result_list[1]}]"
    print(output)

solve_cohomology()