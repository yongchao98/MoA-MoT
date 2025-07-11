def solve_scale_cardinalities():
    """
    This function calculates the cardinalities of the specified group quotients.

    The problem is solved as follows:
    1.  The initial object A in the category of scales is identified as the identity map on integers, id: Z -> Z. Its associated prescale is Z.
    2.  The terminal object B is identified as the zero map Z -> {0}, which is the terminal object in the larger category containing scales. Its prescale is the trivial group {0}.
    3.  The scale S is the inclusion of Z into the hyperreals, *R. Its prescale is *R.
    4.  The quotients are formed from the prescales and the images of the unique canonical maps between them.
        - S/A is *R / Z.
        - B/S is {0} / {0}.
        - B/A is {0} / {0}.
    5.  The cardinalities are computed:
        - |S/A| = |*R/Z|. The cardinality of the hyperreals |*R| is 2^aleph_0 = Beth_1. By Lagrange's theorem for cardinals, |*R/Z| = Beth_1 / aleph_0 = Beth_1.
        - |B/S| = |{0}/{0}| = 1.
        - For H_1(B/A, Q), the space B/A is {0}/{0}, which is a single point. The first homology group of a point, H_1(pt, Q), is the trivial group {0}, which has cardinality 1.
    """
    
    # The calculated values are:
    # Cardinality of S/A: Beth_1
    # Cardinality of B/S: 1
    # Cardinality of H_1(B/A, Q): 1
    
    ans1 = "Beth_1"
    ans2 = "1"
    ans3 = "1"
    
    print(f"{ans1} {ans2} {ans3}")

solve_scale_cardinalities()