import sympy
from fractions import Fraction

# This computation is intensive and may take a moment.

def solve_integral():
    """
    Calculates the definite integral I = 2**50 * integral from 0 to pi of (cos(3x)sin(x))**50 dx
    using symbolic and numerical methods, and returns the result as a fraction.
    """
    x = sympy.Symbol('x')
    
    # The integrand after simplification
    # Since the power is 50 (even), the absolute value is implicitly handled.
    integrand = (2 * sympy.sin(x) * sympy.cos(3*x))**50
    
    # The integral to be calculated
    integral_val = sympy.integrate(integrand, (x, 0, sympy.pi))
    
    # The result from sympy will be in terms of pi.
    # We want to extract the rational coefficient.
    rational_part = integral_val / sympy.pi
    
    # Convert the sympy expression for the rational part to a standard fraction
    # N() evaluates the expression to a floating point number.
    # We use a high precision to ensure accuracy for the fraction conversion.
    float_val = rational_part.evalf(100)
    
    # Convert the float to a fraction. The limit_denominator is important for finding
    # a nice fractional representation of a high-precision float.
    final_fraction = Fraction(float_val).limit_denominator(10**50)
    
    # As requested by the prompt "output each number in the final equation!"
    # The final equation is I = (numerator/denominator) * pi
    print(f"The simplified integral is I = 2^50 * integral from 0 to pi of (cos(3*x)*sin(x))^50 dx.")
    print(f"The computed value of the integral is ({final_fraction.numerator}/{final_fraction.denominator}) * pi.")
    print(f"The fractional part is:")
    print(f"{final_fraction.numerator}/{final_fraction.denominator}")

solve_integral()