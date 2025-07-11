import numpy as np

def solve_frequency_correction():
    """
    This function calculates the coefficient of the 3rd term of the nonlinear
    frequency correction for the Rayleigh-Plesset equation.

    The method of multiple scales is applied to the Rayleigh-Plesset equation,
    perturbed as R = 1 + epsilon * x. The second-order frequency correction omega_2
    is found to be proportional to a polynomial in gamma, which we call C2.

    omega_2 is proportional to -a^2 / (2*omega_0) * C2(gamma)

    C2(gamma) is the sum of three main contributions arising from the
    O(epsilon^2) secularity condition:
    1. TA: From the quadratic nonlinearities in the original equation (x*x_ddot, x_dot^2, x^2).
    2. TB: From the interaction between the fundamental mode (x0) and the first
            harmonic (x1) in the x^2 term.
    3. TC: From the cubic nonlinearity (x^3) in the original equation.

    This script calculates the polynomial C2 by summing these three parts and
    then extracts the coefficient of the third term.
    """

    # Define gamma as a polynomial variable for symbolic-like calculations
    gamma = np.poly1d([1, 0], variable='gamma')

    # Linear oscillation frequency squared
    omega0_sq = 3 * gamma

    # Coefficients from the solution for the first harmonic x1 = K2*A^2*e^(i2*w0*T0) + K0*A*A_bar + c.c.
    # where A is the complex amplitude of the fundamental oscillation.
    K0 = 3 * gamma
    K2 = -(1 + 0.5 * gamma)

    # 1. Calculate the contribution from the quadratic terms (TA)
    # This term combines contributions from x*x_ddot, x_dot^2, and their interactions with x1
    # The combined coefficient for the secular term is omega0_sq * (K0 - K2)
    TA = omega0_sq * (K0 - K2)

    # 2. Calculate the contribution from the x^2 nonlinearity interaction (TB)
    # The coefficient of the secular term is 3*gamma*(3*gamma+1)*(K0+K2)
    TB_prefix = 3 * gamma * (3 * gamma + 1)
    TB = TB_prefix * (K0 + K2)

    # 3. Calculate the contribution from the cubic term x^3 (TC)
    # The coefficient of the secular term is -3*gamma*(3*gamma+1)*(3*gamma+2)/2
    TC = -1.5 * gamma * (3 * gamma + 1) * (3 * gamma + 2)

    # The full coefficient C2 is the sum of these parts
    C2 = TA + TB + TC

    # Output the results
    print("The calculation of the nonlinear frequency correction involves a polynomial in gamma, C2.")
    print("The final polynomial C2(gamma) is the sum of three contributions TA, TB, and TC.")
    print("-" * 50)
    print(f"Contribution from quadratic terms, TA(gamma) = {TA}")
    print(f"Contribution from x^2 interaction, TB(gamma) = {TB}")
    print(f"Contribution from cubic term, TC(gamma) = {TC}")
    print("-" * 50)

    # The final equation for C2
    c3, c2, c1 = C2.coeffs
    print("The full polynomial C2 is expressed as:")
    print(f"C2(gamma) = ({c3}) * gamma^3 + ({c2}) * gamma^2 + ({c1}) * gamma")
    print("-" * 50)

    # The question asks for the 3rd term of the correction. We interpret this as
    # the coefficient of the third term in the polynomial C2, ordered by powers of gamma.
    third_term_coeff = c1
    
    print("The terms of the polynomial C2 are:")
    print(f"1st term: ({c3}) * gamma^3")
    print(f"2nd term: ({c2}) * gamma^2")
    print(f"3rd term: ({c1}) * gamma")
    print("\nThe 3rd term of the nonlinear correction coefficient is therefore the coefficient of gamma.")
    print(f"\nCalculated value of the 3rd term's coefficient: {third_term_coeff}")

if __name__ == '__main__':
    solve_frequency_correction()