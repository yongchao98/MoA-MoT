def solve_lattice_questions():
    """
    Solves the three lattice theory questions based on mathematical theorems and properties.
    """

    # Part (a): Is it true that an even unimodular lattice of rank 12 can have farness exactly 2?
    # A fundamental theorem states that the rank 'n' of an even unimodular lattice must be a multiple of 8.
    rank_a = 12
    # We check if 12 is divisible by 8.
    if rank_a % 8 == 0:
        # If it were possible, we would need further analysis.
        # But since it's not, no such lattice exists.
        answer_a = "Yes"
    else:
        # If no such lattice exists, the statement is false.
        answer_a = "No"

    # Part (b): Suppose L is an odd unimodular lattice of rank 14 with far(L) = 3.
    # Can L have a vector x such that x.x = 0 (mod 6) and x is a 3-primitive vector?
    # This is a deep theoretical question. The existence of such lattices is known.
    # Local-global principles for representation of integers by lattices of rank >= 5
    # suggest that since such a vector can exist locally (in p-adic completions of the lattice),
    # it can also exist globally (in the lattice itself). There's no known theorem forbidding it.
    answer_b = "yes"

    # Part (c): If an even unimodular lattice L in R^24 has a visible root system of type D_24,
    # what is the smallest d for which L can be a d-neighbor of Z^24?
    # The lattice L with root system D_24 is a specific Niemeier lattice,
    # which can be constructed as a 2-neighbor of Z^24 via the D_24 lattice.
    # Therefore, its farness 'd' is at most 2.
    d_is_at_most = 2
    # Farness d=1 means L is isometric to Z^24.
    # L is an even lattice, while Z^24 is an odd lattice, so they are not isometric.
    # Thus, d cannot be 1.
    d_cannot_be = 1
    # The smallest possible value for d is therefore 2.
    answer_c = d_is_at_most

    # Print the final answer in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_lattice_questions()