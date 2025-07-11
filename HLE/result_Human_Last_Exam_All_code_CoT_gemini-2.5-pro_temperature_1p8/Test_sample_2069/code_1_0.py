def generate_milp_constraints():
    """
    This function formulates and prints the two additional inequalities required
    to create an exact MILP encoding for the given function f(x).
    """
    
    # M is a sufficiently large positive constant (Big-M).
    # Its value should be chosen to be larger than any expected value of |y-x| or |y|.
    # A common practice is to choose it based on the bounds l and u, but for a general-purpose
    # representation, a sufficiently large number like 1000 is used.
    M = 1000
    
    # The first new inequality.
    # This inequality is derived from the disjunction on the lower bound of y.
    # It sets the lower bound for y to 0 when b=1 (i.e., x >= 0)
    # and makes the constraint non-binding otherwise.
    # Form: y >= 0 - M*(1-b)  =>  y >= M*b - M
    # It also has the critical role of making wrong choices of (a, b) infeasible.
    inequality1 = f"y >= {M} * b - {M}"
    
    # The second new inequality.
    # This inequality sets the lower bound for y to x when b=0 (i.e., x < 0)
    # and makes the constraint non-binding otherwise.
    # Form: y >= x - M*b
    # This complements the first inequality to correctly model y=min(0,x) when a=0
    # and makes other wrong choices of (a, b) infeasible.
    inequality2 = f"y >= x - {M} * b"
    
    print("The two additional inequalities are:")
    print(inequality1)
    print(inequality2)

generate_milp_constraints()