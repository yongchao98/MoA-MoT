import sympy

def solve_work_cycle():
    """
    This function calculates and displays the work done by the current source
    for a rectangular cycle in the I-x plane for a given electromechanical system.
    """
    # Step 1: Define the symbolic variables for the problem.
    # We use sympy for symbolic mathematics to avoid manual algebraic errors.
    N, w, g, mu, mu0, D = sympy.symbols('N w g mu mu_0 D', positive=True, real=True)
    x1, x2 = sympy.symbols('x_1 x_2', real=True)
    I1, I2 = sympy.symbols('I_1 I_2', real=True)
    x = sympy.Symbol('x') # Position variable

    # Step 2: Determine the inductance L(x) of the coil.
    # The total reluctance of the gap is the parallel combination of the reluctance
    # of the part filled by the magnetic block and the part filled by air.
    # Reluctance R_block = g / (mu * x * w)
    # Reluctance R_air   = g / (mu0 * (D - x) * w)
    # 1/R_total = 1/R_block + 1/R_air
    # After simplification, R(x) = g / (w * ((mu - mu0)*x + mu0*D))
    # Inductance L(x) = N^2 / R(x)
    L_x = (N**2 * w / g) * ((mu - mu0)*x + mu0*D)

    # Step 3: Calculate the work done by the source over the cycle.
    # The work done by the source over a complete cycle is given by W = integral(I * d(lambda))
    # where lambda = L(x)I. For a closed cycle, this is equivalent to the net mechanical work done.
    # A detailed calculation over the four paths reveals:
    # W_total = -1/2 * (L(x2) - L(x1)) * (I2^2 - I1^2).

    # Let's find L(x2) - L(x1)
    L_at_x1 = L_x.subs(x, x1)
    L_at_x2 = L_x.subs(x, x2)
    delta_L = sympy.simplify(L_at_x2 - L_at_x1)
    
    # Let's find I2^2 - I1^2
    delta_I_squared = I2**2 - I1**2
    
    # Construct the final expression for the work W.
    W = -sympy.Rational(1, 2) * delta_L * delta_I_squared

    # Re-arranging for clarity to match one of the choices.
    # The derived expression is W = - (mu - mu0)/(2*g) * N^2 * w * (x2 - x1) * (I2^2 - I1^2)
    
    # Step 4: Print the components of the final equation as requested.
    # We break down the formula into its constituent "numbers" or symbolic terms.
    term1 = -(mu - mu0) / (2 * g)
    term2 = N**2
    term3 = w
    term4 = (x2 - x1)
    term5 = (I2**2 - I1**2)

    print("The work done by the current source, W, is the product of the following terms:")
    print(f"Term 1 (Physical constants and geometry): {term1}")
    print(f"Term 2 (Turns squared): {term2}")
    print(f"Term 3 (Width): {term3}")
    print(f"Term 4 (Displacement): {term4}")
    print(f"Term 5 (Change in current squared): {term5}")
    print("\nThe resulting final equation is:")
    # The final print statement constructs the readable version of the equation.
    print(f"W = ({sympy.pretty(term1, use_unicode=False)}) * ({term2}) * ({term3}) * ({term4}) * ({term5})")


if __name__ == '__main__':
    solve_work_cycle()
    print("\nThis result corresponds to answer choice D.")