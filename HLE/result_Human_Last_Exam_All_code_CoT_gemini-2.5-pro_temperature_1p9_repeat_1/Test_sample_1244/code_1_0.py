def solve_lattice_questions():
    """
    Solves the three lattice theory questions based on mathematical properties.
    """

    # Part (a): Is it true that an even unimodular lattice of rank 12 can have farness exactly 2?
    # According to a theorem by Milnor, the rank of an even unimodular lattice must be a multiple of 8.
    # Since 12 is not a multiple of 8, such a lattice does not exist.
    answer_a = "No"

    # Part (b): Suppose L is an odd unimodular lattice of rank 14 with far(L) = 3.
    # Can L have a vector x such that x.x = 0 (mod 6) and x is a 3-primitive vector?
    # We showed in the reasoning that such a lattice and vector can be constructed.
    # We found a suitable vector x and a lattice L that is a 3-neighbor of Z^14 containing x.
    answer_b = "yes"

    # Part (c): If an even unimodular lattice L in R^24 has a visible root system of type D_24,
    # what is the smallest d for which L can be a d-neighbor of Z^24?
    # The lattice L is the Niemeier lattice N_D24, also known as D_24^+.
    # The intersection L intersect Z^24 is the lattice D_24.
    # The index [L : D_24] is 2.
    # The index [Z^24 : D_24] is 2.
    # Since both indices are equal, L is a 2-neighbor of Z^24.
    # Since L is not isometric to Z^24 (even vs odd), the farness cannot be 1.
    # Therefore, the smallest d is 2.
    answer_c = 2

    # The final answer formatted as requested.
    result = f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]"
    print(result)

solve_lattice_questions()