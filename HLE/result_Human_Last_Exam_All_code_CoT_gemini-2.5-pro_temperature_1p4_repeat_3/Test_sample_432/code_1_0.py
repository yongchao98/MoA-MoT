def solve_cardinality_problem():
    """
    This function identifies which of a given list of infinite sets have the same
    cardinality as the interval [0, 1] and prints the answer.
    """

    # The cardinality of [0, 1] is the cardinality of the continuum, 'c'.
    # We will identify all sets from the list with cardinality 'c'.
    # The reasoning for each set's cardinality is provided in the comments.
    # The numbers in the 'equations' of cardinal arithmetic (e.g., the '2' in R^2)
    # are mentioned in the comments.

    cardinalities = {
        # --- Sets with cardinality 'c' (the continuum) ---
        'A': 'c',  # |(0, 1)| = c. A simple bijection maps (0,1) to R.
        'D': 'c',  # |R| is the definition of the continuum, c.
        'E': 'c',  # |R \ Q| (irrationals). |R| = |Q|+|R\Q| => c = aleph_0 + |R\Q| => |R\Q|=c.
        'F': 'c',  # |C| = |R^2|. The power is 2. |R^2|=c^2=c.
        'G': 'c',  # |H| = |R^4|. The power is 4. |R^4|=c^4=c.
        'H': 'c',  # |{x: c'(x) = 0}| = |[0,1]\Cantor set| = aleph_0 * c = c.
        'J': 'c',  # |R^N| (infinite sequences of reals). |R^N|=c^aleph_0=c.
        'K': 'c',  # |Z^N| (infinite sequences of integers). |Z^N|=aleph_0^aleph_0=c.
        'M': 'c',  # |R x R| = |R^2|. The power is 2. |R^2|=c^2=c.
        'N': 'c',  # |2^N| (power set of naturals). |2^N|=2^aleph_0=c.
        'O': 'c',  # |2^Q| (power set of rationals). |2^Q|=2^aleph_0=c.

        # --- Sets with countable cardinality 'aleph_0' ---
        'B': 'aleph_0', # |N| = aleph_0.
        'C': 'aleph_0', # |Q| = aleph_0.
        'I': 'aleph_0', # Set of all finite strings over a finite alphabet is countable.
        'L': 'aleph_0', # |N x N| = aleph_0^2 = aleph_0. The power is 2.

        # --- Sets with cardinality '2^c' (greater than continuum) ---
        'P': '2^c', # |2^C| (power set of complex numbers) = 2^c.
        'Q': '2^c'  # |2^(2^Z)| = 2^(2^aleph_0) = 2^c.
    }

    # Filter for sets with the target cardinality 'c'
    target_cardinality = 'c'
    correct_options = [option for option, card in cardinalities.items() if card == target_cardinality]

    # Sort the options alphabetically
    correct_options.sort()

    # Join into a single string and print
    final_answer = "".join(correct_options)
    print(final_answer)

solve_cardinality_problem()