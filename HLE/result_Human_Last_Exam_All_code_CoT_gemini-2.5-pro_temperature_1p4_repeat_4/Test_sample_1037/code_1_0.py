import math

def solve():
    """
    This function calculates the cardinalities of S/A, B/S, and H_1(B/A, Q).

    1.  |S/A|: This is the cardinality of the quotient group of the hyperreals by the integers, |*R / Z|.
        The cardinality of the hyperreals *R is beth_1 (the power of the continuum).
        The cardinality of the integers Z is beth_0.
        The cardinality of the quotient |*R / Z| is beth_1.
        We represent this as "Beth_1".

    2.  |B/S|: This is the cardinality of the quotient group G_B / Im(h_SB).
        Based on the analysis of the terminal object B, its underlying group G_B is the trivial group {0}.
        The canonical map h_SB from S to B must be the zero map. Its image is {0}.
        The quotient is {0}/{0}, which is the trivial group {0}.
        Its cardinality is 1.

    3.  |H_1(B/A, Q)|: This is the cardinality of the first homology group of B/A.
        The group B/A is G_B / Im(h_AB), which is {0}/{0} = {0}.
        The first homology group of the trivial group {0} with any coefficients is the trivial group.
        So H_1({0}, Q) = {0}.
        Its cardinality is 1.

    The final answer combines these three parts.
    """

    # The cardinalities determined by the reasoning above.
    card_S_div_A = "Beth_1"
    card_B_div_S = 1
    card_H1_B_div_A = 1
    
    # The problem asks to output the equation, so we print the values.
    # The values are determined by mathematical deduction, not direct computation.
    # This code serves to format and print the final result.
    
    print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")

solve()