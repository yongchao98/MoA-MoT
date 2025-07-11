import sympy

def display_minimax_risk_formula():
    """
    This function calculates and displays the symbolic minimax risk for estimating
    the parameter theta of a Binomial distribution under the specified conditions.

    The derivation shows that for a sufficient statistic S ~ Bin(N, theta),
    the minimax risk is 1 / (4 * (sqrt(N) + 1)**2).
    In this problem, N = n*n.
    """
    # Define n as a symbolic variable
    n = sympy.Symbol('n', positive=True, integer=True)

    # Total number of trials is N = n^2
    N = n**2

    # The general formula for minimax risk for Bin(N, theta)
    minimax_risk_expr = 1 / (4 * (sympy.sqrt(N) + 1)**2)

    # Simplify the expression
    minimax_risk_simplified = sympy.simplify(minimax_risk_expr)

    # Print the final formula in a human-readable format
    print("Based on the derivation, the minimax risk for estimating theta is given by the formula:")
    sympy.pprint(minimax_risk_simplified, use_unicode=True)

    # As requested, output each number in the final equation: 1 / (4*(n+1)**2)
    print("\nThe integer constants in the final equation are:")
    print("Numerator: 1")
    print("Denominator factor: 4")
    print("Term added to n inside parenthesis: 1")
    print("Exponent on the parenthesis: 2")

if __name__ == '__main__':
    display_minimax_risk_formula()
