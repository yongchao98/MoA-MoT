def solve_lattice_problems():
    """
    Solves the three lattice theory problems and prints the answer.
    """

    # Part (a): Is it true that an even unimodular lattice of rank 12 can have farness exactly 2?
    # An even unimodular lattice L of rank 12 must have a signature (n_+, n_-)
    # such that n_+ + n_- = 12 and n_+ - n_- is a multiple of 8.
    # The possible signatures are (10, 2), (6, 6), and (2, 10).
    # The maximal dimension of a positive-definite subspace is n_+.
    # farness(L) = 2 implies L contains a sublattice isometric to 2*Z^12.
    # This sublattice is a positive-definite subspace of dimension 12.
    # However, for a rank 12 even unimodular lattice, the maximum possible n_+ is 10.
    # Since 12 > 10, no such sublattice can exist.
    answer_a = "No"

    # Part (b): Suppose L is an odd unimodular lattice of rank 14 with far(L) = 3.
    # Can L have a vector x such that x.x is congruent to 0 (mod 6) and x is a 3-primitive vector?
    # far(L) = 3 implies L contains a sublattice M isometric to 3*Z^14.
    # L can be constructed from M using a self-dual code C over F_3 of length 14.
    # A vector x in L can be a "glue vector" whose norm x.x is the weight of a codeword in C.
    # The Pless code P_14 is a suitable code, and it contains codewords of weight 6.
    # We can choose a vector x corresponding to such a codeword, so x.x = 6.
    # This satisfies x.x = 6, which is congruent to 0 (mod 6).
    # A vector is 3-primitive if x/3 is not in L.
    # If x is a glue vector with norm 6, x/3 cannot be in L, because that would imply
    # the corresponding codeword is the zero vector, but its weight is 6.
    # Thus, such a vector x can exist.
    answer_b = "yes"

    # Part (c): If an even unimodular lattice L in R^24 has a visible root system of type D_24,
    # what is the smallest d for which L can be a d-neighbor of Z^24?
    # The lattice L is the Niemeier lattice N(D_24).
    # It is constructed as L = D_24 U (h + D_24), where D_24 is the lattice of integer vectors
    # with an even sum of coordinates, and h = (1/2, ..., 1/2).
    # L is a d-neighbor of Z^24 if the intersection M = L intersect Z^24 has index d
    # in both L and Z^24.
    # The intersection is M = D_24.
    # The index of D_24 in Z^24 is 2, since D_24 is the kernel of the sum-of-coordinates-mod-2 map.
    # The index of D_24 in L is also 2, by the construction of L.
    # So, L is a 2-neighbor of Z^24.
    # d=1 is not possible because L is an even lattice and Z^24 is odd, so they are not isometric.
    # Therefore, the smallest d is 2.
    answer_c = 2

    # Formatting the final answer
    final_answer = f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]."
    print(final_answer)

solve_lattice_problems()