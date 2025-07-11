import sympy

def solve_limit():
    """
    This function calculates and displays the symbolic limit based on the derivation.
    """
    # Define k as a symbolic variable, representing the integer k >= 2.
    k = sympy.Symbol('k', integer=True, positive=True)

    # The limit was derived to be k / (k + 1).
    limit_expression = k / (k + 1)

    # The problem asks to output the equation and the numbers in it.
    # The equation is L = k / (k + 1), which can be written as L = (1*k) / (1*k + 1).
    # The numbers in the final equation are the coefficients of k and the constants.
    
    print("The computed limit is a symbolic expression in terms of k.")
    
    # We construct the equation string, showing the numbers (coefficients) explicitly.
    numerator_str = f"1*k"
    denominator_str = f"1*k + 1"
    
    final_equation_str = f"L = ({numerator_str}) / ({denominator_str})"
    
    print(f"The final equation for the limit L is: {final_equation_str}")
    
    # As requested, we output each number in the final equation.
    print("\nThe numbers appearing in the final equation are:")
    print(f"Coefficient of k in the numerator: 1")
    print(f"Coefficient of k in the denominator: 1")
    print(f"Constant term in the denominator: 1")

if __name__ == '__main__':
    solve_limit()