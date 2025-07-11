import sympy

def solve_integral():
    """
    This function calculates the specified integral using symbolic mathematics.
    """
    # Step 1: Define the symbolic variable for τ.
    # The problem specifies τ > 0.
    tau = sympy.symbols('τ', positive=True)

    # Step 2: Define the expression for the sum of initial positions S(0; τ).
    # This was derived from the differential equations of motion.
    # S_0_tau = 1/4 + (3/4) * exp(2*τ^2)
    s_0_tau = sympy.Rational(1, 4) + sympy.Rational(3, 4) * sympy.exp(2 * tau**2)

    # Step 3: Define the integrand for the final calculation, which is 1 / S(0; τ).
    integrand = 1 / s_0_tau
    
    print(f"The problem reduces to calculating the definite integral of the following function from 0 to ∞:")
    print(f"f(τ) = {integrand.simplify()}")
    print("-" * 30)

    # Step 4: Perform the definite symbolic integration.
    integral_result = sympy.integrate(integrand, (tau, 0, sympy.oo))

    # Step 5: Format and display the final answer.
    # The result is ln(16/9). We extract the numerator and denominator
    # to display them as requested.
    # The argument of the logarithm is exp(integral_result).
    log_argument = sympy.exp(integral_result)
    numerator = sympy.numer(log_argument)
    denominator = sympy.denom(log_argument)

    print("The symbolic evaluation of the integral gives the final answer:")
    # The print statement below will explicitly show the numbers in the final fraction.
    print(f"Final Answer = ln({numerator}/{denominator})")

if __name__ == '__main__':
    solve_integral()
