def solve_dispersion_point_problem():
    """
    This function solves the mathematical problem about the maximum cardinality
    of the set of dispersion points in a compact connected metric space.
    
    The reasoning is as follows:
    1.  It can be proven by contradiction that a compact connected metric space
        can have at most one dispersion point. The proof relies on the properties
        of irreducible continua. If we assume two dispersion points d1 and d2 exist,
        we can find an irreducible continuum K between them. K\{d1} must be connected
        (a property of irreducible continua) but also totally disconnected (as a subset
        of X\{d1}). This forces K to be the two-point set {d1, d2}, which is not
        connected in a metric space, leading to a contradiction.
    2.  Examples exist, such as the Knaster-Kuratowski fan, which are compact
        connected metric spaces with exactly one dispersion point.
    
    Combining these two facts, the maximum possible cardinality is 1.
    """
    
    # The maximum cardinality of the set of dispersion points.
    maximum_cardinality = 1
    
    # The final equation is simply stating this result.
    # The problem asks to output each number in the final equation.
    # We will print the equation: "Maximum Cardinality = 1"
    
    print(f"Maximum Cardinality = {maximum_cardinality}")

solve_dispersion_point_problem()