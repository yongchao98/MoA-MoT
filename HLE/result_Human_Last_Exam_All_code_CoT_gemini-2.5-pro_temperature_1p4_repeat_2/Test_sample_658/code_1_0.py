import math

def get_asymptotic_formula():
    """
    This function returns the asymptotic formula for d_{B,delta}.
    The final expression is built using string formatting.
    """
    B = "B"
    L = "L"
    
    # The formula is a combination of the two asymptotic regimes.
    # We use 'sqrt' for square root and 'log' for natural logarithm.
    # The constant '2' in log(L/B + 2) is to ensure the argument is always > 1,
    # as L/B can be close to 0 but is always positive.
    
    term1 = f"sqrt({B} * {L})"
    term2_numerator = f"{L}"
    term2_denominator = f"log({L} / {B} + 2)"
    
    formula = f"{term1} + {term2_numerator} / ({term2_denominator})"
    
    # We are asked to output each number, which is done here with the '2'.
    print(formula)

get_asymptotic_formula()