def solve_milp_encoding():
    """
    This function formulates and prints the two additional inequalities required
    to create an exact MILP encoding of the given function f(x).
    """

    # The original constraints are given for variables x, y, a, with l <= x <= u.
    # a is a binary variable (0 or 1).
    #
    # The function f(x) is:
    # f(x) = x - 1, if x >= 1
    # f(x) = 0,     if 0 <= x < 1
    # f(x) = x,     if x < 0
    #
    # We introduce a new binary variable 'b' to model the piecewise nature of f(x) for x < 1.
    # The logic is as follows:
    # If a=0, b=1 (intended for 0 <= x < 1): we need to enforce y = 0.
    # If a=0, b=0 (intended for x < 0): we need to enforce y = x.
    #
    # The existing constraints for a=0 give: x-1 <= y <= min(0,x).
    # We need to add constraints to select the correct value for y based on b.
    #
    # New Inequality 1: This inequality helps enforce y >= 0 when b=1.
    # Using the lower bound 'l' for x as a parameter for a tight formulation.
    # The inequality is y >= l - l*b.
    # if b=1: y >= l - l*1 => y >= 0.
    # if b=0: y >= l - l*0 => y >= l (a valid, non-restrictive bound).
    inequality1 = "y >= l - l*b"

    # New Inequality 2: This inequality helps enforce y >= x when b=0.
    # Using the upper bound 'u' for x as a parameter for a tight formulation.
    # The inequality is y >= x - u*b.
    # if b=0: y >= x - u*0 => y >= x.
    # if b=1: y >= x - u*1 => y >= x - u (a valid, non-restrictive bound, as y's known lower bound is l).
    inequality2 = "y >= x - u*b"

    # The problem asks to format the answer as a comma-separated list of inequalities.
    # The terms are already distributed (no brackets).
    print(f"{inequality1}, {inequality2}")

solve_milp_encoding()