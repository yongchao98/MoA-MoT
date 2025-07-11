import sympy
import numpy as np

def calculate_limit():
    """
    This function calculates the limit of p_n as n approaches infinity by modeling
    the hitting probability using the walk's asymptotic angular distribution.
    """
    # Define symbolic variables for the calculation
    n = sympy.Symbol('n', positive=True, real=True)
    c = sympy.Symbol('c', positive=True, real=True) # Represents a constant of order 1
    theta = sympy.Symbol('theta', real=True)

    # The probability density function (PDF) for the limiting angle theta is (sin(theta)^2)/pi.
    # For small theta, sin(theta) is approximately theta.
    pdf_approx = theta**2 / sympy.pi

    # The walk hits the target if its limiting angle is within a small interval [-delta, delta]
    # around 0. The width of this interval, delta, is of order 1/sqrt(n).
    delta = c / sympy.sqrt(n)

    # We approximate p_n by integrating the PDF over this interval.
    p_n_approx = sympy.integrate(pdf_approx, (theta, -delta, delta))
    
    # Calculate the limit of this expression as n goes to infinity.
    limit_p_n = sympy.limit(p_n_approx, n, sympy.oo)

    print("Step-by-step calculation:")
    print("1. The probability density of the limiting angle theta is f(theta) = (sin(theta)^2)/pi.")
    print("2. For the walk to hit the target, theta must be close to 0.")
    print("3. We approximate p_n by the integral of f(theta) over an interval [-c/sqrt(n), c/sqrt(n)].")
    print(f"   p_n ≈ Integral from {-delta} to {delta} of ({pdf_approx}) d_theta")
    
    # Let's show the parts of the resulting equation
    numerator, denominator = p_n_approx.as_numer_denom()
    
    # To satisfy the output format requirement of showing each number
    num_parts = numerator.as_ordered_terms()
    den_parts = denominator.as_ordered_terms()
    
    print(f"4. The integral evaluates to: ({num_parts[0]}) / ({den_parts[0]} * {den_parts[1]})")
    
    num_coeff = num_parts[0] / c**3
    den_coeff_1 = den_parts[0] / sympy.pi
    den_coeff_2 = den_parts[1] / n**sympy.S('3/2')
    
    print(f"   The equation for p_n is approximately ({int(num_coeff)}*c^3) / ({int(den_coeff_1)}*pi*n^({str(den_coeff_2*2)}))")

    print(f"5. The limit of this expression as n -> infinity is:")
    print(limit_p_n)

if __name__ == '__main__':
    calculate_limit()
