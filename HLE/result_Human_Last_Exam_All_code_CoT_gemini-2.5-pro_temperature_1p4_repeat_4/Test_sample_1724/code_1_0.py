import sympy

def solve_frequency_correction():
    """
    Calculates the polynomial P(gamma) that determines the second-order
    frequency correction w2 for the Rayleigh-Plesset equation.

    The secularity condition at order epsilon^3 leads to the expression:
    2*w0*w2 = (3*gamma / 8) * P(gamma)

    This function calculates the coefficients of P(gamma) = A*gamma^2 + B*gamma + C
    and identifies the third term, C.
    """

    # Define gamma as a symbol
    gamma = sympy.symbols('gamma')

    # The polynomial P(gamma) is a sum of four parts, arising from different
    # terms in the perturbation expansion at order epsilon^3. We define them here.
    # These terms are derived from a detailed Poincaré-Lindstedt analysis.

    # Part 1: from omega_0^2 * (R1*R2'' + R2*R1'')
    p1 = 10 - gamma
    
    # Part 2: from 3 * omega_0^2 * (R1'*R2')
    p2 = -6 * (gamma + 2)

    # Part 3: from -3*gamma*(3*gamma+1) * (R1*R2)
    p3 = -(3 * gamma + 1) * (5 * gamma - 2)
    
    # Part 4: from (gamma*(3*gamma+1)*(3*gamma+2)/2) * R1^3
    p4 = (3 * gamma + 1) * (3 * gamma + 2)
    
    # The full polynomial P(gamma)
    P_gamma = sympy.expand(p1 + p2 + p3 + p4)

    # Extract the coefficients of the polynomial P(gamma) = A*gamma^2 + B*gamma + C
    coeff_A = P_gamma.coeff(gamma, 2)
    coeff_B = P_gamma.coeff(gamma, 1)
    coeff_C = P_gamma.coeff(gamma, 0)
    
    print(f"The determining polynomial for the frequency correction is P(gamma) = {P_gamma}")
    print("\nThis polynomial can be written as A*gamma^2 + B*gamma + C")
    print("The final equation for the terms of the polynomial is:")
    print(f"({coeff_A}) * gamma**2 + ({coeff_B}) * gamma + ({coeff_C})")

    # The "3rd term" is interpreted as the constant term, C.
    third_term = coeff_C
    
    print(f"\nThe third term of this polynomial is the constant term, C = {third_term}")
    print("\nFinal Answer:")
    print(f"<<<{third_term}>>>")


if __name__ == '__main__':
    solve_frequency_correction()
