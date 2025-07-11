def solve_lattice_problems():
    """
    This function provides the solution to the three-part lattice theory problem.
    The reasoning for each part is commented below.
    """

    # Part (a): Is it true that an even unimodular lattice of rank 12 can have farness exactly 2?
    # A positive-definite even unimodular lattice must have a rank that is a multiple of 8.
    # Since n=12 is not a multiple of 8, such a lattice does not exist.
    # A lattice must be positive-definite to be a neighbor of the standard lattice Z^n.
    answer_a = "No"

    # Part (b): Suppose L is an odd unimodular lattice of rank 14 with far(L) = 3.
    # Can L have a vector x such that x.x is congruent to 0 (mod 6) and x is a 3-primitive vector?
    # For a unimodular lattice L, x being 3-primitive means y = x/3 is a primitive vector in L.
    # The condition x.x = 0 (mod 6) implies 9*(y.y) = 0 (mod 6), which requires y.y to be even.
    # An odd lattice can contain primitive vectors with an even norm.
    # A construction of such a lattice L confirms the existence of a primitive vector y with y.y=2.
    answer_b = "yes"

    # Part (c): If an even unimodular lattice L in R^24 has a visible root system of type D_24,
    # what is the smallest d for which L can be a d-neighbor of Z^24?
    # The specified lattice L is the Niemeier lattice with root system D_24, which is D_24^+.
    # L is a d-neighbor of Z^24 if K = L intersect Z^24 has index d in both lattices.
    # The intersection is K = D_24.
    # The index d = [Z^24 : D_24] is 2.
    # Since L is not isometric to Z^24, d must be greater than 1.
    # Therefore, the smallest d is 2.
    answer_c = 2

    # Print the final answer in the required format.
    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]")

solve_lattice_problems()