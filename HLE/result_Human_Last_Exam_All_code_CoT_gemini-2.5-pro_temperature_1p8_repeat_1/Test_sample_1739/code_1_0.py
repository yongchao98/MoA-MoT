import sympy

def solve_rayleigh_plesset_correction():
    """
    This script calculates the 3rd term of the nonlinear frequency correction for the Rayleigh-Plesset equation.

    The nonlinear frequency correction is expressed as a coefficient C1 in the relation:
    omega^2 = omega_0^2 * (1 - C1 * a^2)
    where 'a' is the dimensionless amplitude of the oscillation.

    The script calculates C1, which is a quadratic polynomial in the polytropic index gamma:
    C1 = A*gamma^2 + B*gamma + C

    The "3rd term" is interpreted as the constant term C in this polynomial.
    """

    # Define gamma as a symbolic variable
    gamma = sympy.Symbol('g')

    # Linear oscillation frequency squared
    omega_0_sq = 3 * gamma

    # Coefficients for the first-order solution x1.
    # x1 = K0*|A|^2 + K2*(A^2*exp(2*i*w0*T0) + c.c.), where x0 = A*exp(i*w0*T0) + c.c.
    # These are derived from the O(epsilon) equation's forcing terms.
    K0 = 3 * gamma
    K2 = -(gamma + 2) / 2

    # S is the coefficient of the secular term |A|^2*A*exp(i*w0*T0) in the O(epsilon^2) equation.
    # S = omega_0_sq*(K0-K2) + (9*gamma**2+3*gamma)*(K0+K2) - 3*(9/2*gamma**3+9/2*gamma**2+gamma)
    
    term1 = omega_0_sq * (K0 - K2)
    term2 = (9 * gamma**2 + 3 * gamma) * (K0 + K2)
    term3 = -3 * (sympy.Rational(9, 2) * gamma**3 + sympy.Rational(9, 2) * gamma**2 + gamma)

    S = term1 + term2 + term3

    # The frequency correction coefficient C1 is derived from the secularity condition
    # which relates the frequency shift to S: C1 = S / (4 * omega_0_sq)
    C1 = sympy.expand(S / (4 * omega_0_sq))

    # C1 is a quadratic polynomial in gamma. We extract the coefficients.
    poly_C1 = sympy.Poly(C1, gamma)
    coeffs = poly_C1.all_coeffs()
    
    A = coeffs[0]
    B = coeffs[1]
    C = coeffs[2]

    # Output the components of the equation for C1
    print("The coefficient of the nonlinear frequency correction, C1, is a quadratic polynomial in the polytropic index gamma.")
    print(f"The equation for C1 is: ({A}) * gamma**2 + ({B}) * gamma + ({C})")
    
    # Print the identified 3rd term
    print("\nThe 3rd term of this expression is the constant term C.")
    print(C)


if __name__ == "__main__":
    solve_rayleigh_plesset_correction()