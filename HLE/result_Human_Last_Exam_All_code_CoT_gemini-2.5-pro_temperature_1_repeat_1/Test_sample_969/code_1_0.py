def solve_sequence():
    """
    This function determines the next 4 elements of the sequence based on its underlying pattern.
    The pattern is a run-length encoding where:
    1. The sequence of values (V) bounces between 1 and 3 (3,2,1,2,3,2,1,...).
    2. The sequence of lengths (L) follows a discovered pattern (1,1,1,1,3,3,1,1,1,...).
    The provided sequence is a partial representation, and we complete it based on this pattern.
    """

    # The established pattern for run lengths, discovered by analyzing the completed sequence.
    pattern_l = [1, 1, 1, 1, 3, 3, 1, 1, 1]

    # The given sequence has 9 elements. Let's determine the next 4.
    # s_1 to s_9 are given. We need s_10, s_11, s_12, s_13.

    # From the analysis, the 6th run (of value 2) is incomplete.
    # Its length L_6 should be 3 according to the pattern. We have two 2s.
    s_10 = 2

    # The 7th run has value V_7 = 1 and length L_7 = 1.
    s_11 = 1

    # The 8th run has value V_8 = 2 and length L_8 = 1.
    s_12 = 2

    # The 9th run has value V_9 = 3 and length L_9 = 1.
    s_13 = 3

    next_four_elements = [s_10, s_11, s_12, s_13]

    print(f"The next 4 elements of the sequence are: {next_four_elements[0]} {next_four_elements[1]} {next_four_elements[2]} {next_four_elements[3]}")

solve_sequence()