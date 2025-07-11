import sympy

def display_asymptotic_behavior():
    """
    This function calculates and prints the corrected large-distance behavior
    of the solution omega to the given PDE.
    """

    # Define the symbolic variables used in the equation
    r = sympy.Symbol('r', positive=True)
    theta = sympy.Symbol('theta', real=True)
    A = sympy.Symbol('A')
    B = sympy.Symbol('B')
    omega = sympy.Function('omega')(r, theta)

    # --- Construct the components of the asymptotic formula ---

    # The original solution for A=B=0 behaves like r**(-1/2). This term persists.
    base_power = sympy.Rational(-1, 2)

    # The presence of non-zero A and B adds a corrector to the power of r.
    # This correction is derived from the potential Phi and the WKB analysis.
    power_corrector = A * (1 - sympy.cos(theta)) + B * sympy.sin(theta)

    # The total exponent for the r term in the solution
    total_power = power_corrector + base_power

    # An additional angular phase factor is introduced due to the rotational part
    # of the vector field u, related to the parameter B.
    angular_factor = sympy.exp(B * theta)

    # The dominant exponential decay term is the same as in the A=B=0 case.
    exponential_decay = sympy.exp(-r * (1 - sympy.cos(theta)))

    # Combine all parts to form the final asymptotic expression for omega
    asymptotic_omega = r**total_power * angular_factor * exponential_decay

    # --- Print the final result ---

    print("The corrected large-distance behavior of ω(r, θ) is given by the asymptotic formula:")
    
    # We use sympy's pretty printer for a clear, typeset-style mathematical output.
    # This will clearly show all terms and numbers, including the coefficients 1 and -1/2.
    final_equation = sympy.Eq(omega, asymptotic_omega, evaluate=False)
    sympy.pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    display_asymptotic_behavior()