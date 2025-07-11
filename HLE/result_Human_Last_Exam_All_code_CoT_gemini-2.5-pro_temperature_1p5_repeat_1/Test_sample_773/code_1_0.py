def final_answer_formula():
    """
    This function prints the symbolic formula for the requested total mass.
    The formula involves symbolic variables n, q, q_v and the zeta function Z.
    """

    # The problem asks for a formula for the total mass. The derivation shows
    # the answer depends on n, q, q_v, and the Dedekind zeta function Z of the ring R.
    # The final expression is: (q_v / (q_v - 1)) * Z(2) * Z(3) * ... * Z(n).

    # The prompt requires outputting each number in the final equation.
    # The numbers involved are:
    # 1: in the term (q_v - 1)
    # 2, 3, ..., n: as arguments to the zeta function Z.

    # We construct and print the equation string to display this result clearly.
    equation_string = f"Total Mass = (q_v / (q_v - 1)) * Z(2) * Z(3) * ... * Z(n)"

    print(equation_string)

final_answer_formula()