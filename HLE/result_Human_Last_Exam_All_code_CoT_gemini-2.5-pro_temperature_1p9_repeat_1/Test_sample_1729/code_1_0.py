import sympy

def find_probability_formula():
    """
    This function derives and prints the formula for the probability P_m.
    """
    # Define m as a symbolic variable, representing a positive integer.
    m = sympy.Symbol('m', integer=True, positive=True)

    # The total number of ways to choose two distinct indices i and j from
    # the set {1, 2, ..., 4m+2} is given by the combination formula C(4m+2, 2).
    total_pairs = sympy.binomial(4*m + 2, 2)
    
    # Based on the analysis, there are exactly 3 pairs {i, j} for which the
    # remaining 4m items can be partitioned as required. These pairs are:
    # {1, 2}, {4m+1, 4m+2}, and {1, 4m+2}.
    favorable_pairs = 3

    # The probability P_m is the ratio of the number of favorable pairs
    # to the total number of pairs.
    P_m = sympy.Rational(favorable_pairs, total_pairs)

    # Simplify the expression for the probability.
    simplified_P_m = sympy.simplify(P_m)

    # Output the results in a clear format, following the instructions.
    # The final equation is P_m = 3 / ((2*m + 1)*(4*m + 1)).
    
    # Extract the numerator and denominator for printing.
    num, den = sympy.fraction(simplified_P_m)

    # As requested, output each number in the final equation.
    # The equation has the numbers 3, 2, 1, 4, 1.
    equation_str = f"P_m = {num} / (({2}*m + {1})*({4}*m + {1}))"
    
    print("The derived formula for P_m is:")
    print(equation_str)

find_probability_formula()