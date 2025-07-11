def solve_scale_cardinalities():
    """
    This function determines the cardinalities based on the mathematical derivation.

    1.  |S/A|: The quotient of the hyperreals by the integers.
        The cardinality of the hyperreals |*R| is the continuum, Beth_1.
        Quotienting by the countable set of integers Z does not change the cardinality.
        So, |S/A| = Beth_1.

    2.  |B/S|: The quotient of the terminal object's group by the image of the hyperreals.
        The terminal object B is the trivial scale (Z -> {0}), so its group is G_B = {0}.
        The quotient is {0} / im(h), which is the trivial group {0}.
        The cardinality is 1.

    3.  |H_1(B/A, Q)|: The cardinality of the first homology group of B/A.
        The space B/A is G_B / G_A = {0} / Z, which is the trivial group {0}.
        This space is a single point. The first homology group of a point is 0.
        The cardinality is 1.
    """

    # The derived cardinalities
    card_S_div_A = "Beth_1"
    card_B_div_S = "1"
    card_H1_B_div_A = "1"

    # The problem asks to output each number in the final equation.
    print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")

solve_scale_cardinalities()