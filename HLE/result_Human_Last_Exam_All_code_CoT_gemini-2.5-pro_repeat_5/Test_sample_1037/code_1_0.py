def solve_scale_cardinalities():
    """
    This function calculates the cardinalities requested in the problem.

    The problem asks for three values:
    1. The cardinality of S/A.
    2. The cardinality of B/S.
    3. The cardinality of H_1(B/A, Q).

    Our step-by-step analysis determined these values to be:
    1. |S/A| = Beth_1
    2. |B/S| = 1
    3. |H_1(B/A, Q)| = 1
    """

    # Cardinality of S/A, which is |*R / Z|
    card_S_div_A = "Beth_1"

    # Cardinality of B/S, which is |{0} / {0}|
    card_B_div_S = 1

    # Cardinality of H_1(B/A, Q), which is |H_1({0}, Q)|
    card_H1_B_div_A = 1

    # The problem asks for the output in the format "value1 value2 value3"
    # and to output each number in the final equation.
    print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")

solve_scale_cardinalities()