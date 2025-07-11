def solve_feynman_diagram_count():
    """
    This function determines the number of non-vanishing Feynman diagrams of order 2n=6
    for the electron propagator in quantum electrodynamics.

    The sequence a(n) represents the number of one-particle-irreducible (1PI)
    Feynman diagrams for the electron self-energy at n-loop order.
    The order of the diagram is 2n.
    
    The first few terms of the sequence are:
    a(1) = 1  (for 1-loop, order 2)
    a(2) = 7  (for 2-loops, order 4)
    a(3) = 72 (for 3-loops, order 6)
    """

    # We use a dictionary to store the known values of the sequence a(n).
    a_n_sequence = {
        1: 1,
        2: 7,
        3: 72,
    }

    n = 3
    result = a_n_sequence.get(n)

    # The final output is an equation showing each number, as requested.
    if result is not None:
        print(f"a({n}) = {result}")
    else:
        print(f"The value for a({n}) is not defined in the sequence.")

solve_feynman_diagram_count()