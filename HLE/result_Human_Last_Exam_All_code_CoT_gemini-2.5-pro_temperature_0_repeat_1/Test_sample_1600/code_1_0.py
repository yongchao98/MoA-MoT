def solve_feynman_diagram_count():
    """
    Calculates a(3) from the sequence a(n) for the number of non-vanishing
    Feynman diagrams of order 2n in QED.
    """
    # The sequence a(n) gives the number of non-vanishing Feynman diagrams of
    # order 2n for the electron or photon propagators in quantum electrodynamics.
    # The first few known values of the sequence are stored in a dictionary.
    # The key is 'n' and the value is 'a(n)'.
    a_n_values = {
        1: 1,
        2: 2,
        3: 10,
        4: 74,
        5: 706
    }

    # The problem asks for the value of a(3).
    n = 3

    # Retrieve the value from the dictionary.
    if n in a_n_values:
        result = a_n_values[n]
        # The problem asks to show the numbers in the final equation.
        # We will print the relationship a(n) = result.
        print(f"The number of non-vanishing Feynman diagrams of order 2n is given by the sequence a(n).")
        print(f"For n = {n}, the order is 2*n = {2*n}.")
        print(f"The value is a({n}) = {result}.")
    else:
        print(f"The value for a({n}) is not available in the pre-computed list.")

solve_feynman_diagram_count()