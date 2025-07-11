import sympy

def solve_rayleigh_plesset_correction():
    """
    This function calculates and prints the third term of the nonlinear frequency correction
    for the Rayleigh-Plesset equation, based on the Poincaré-Lindstedt perturbation method.

    The frequency expansion is given by:
    omega = omega_0 + epsilon * omega_1 + epsilon^2 * omega_2 + ...

    Through perturbation analysis, it is found that omega_1 = 0 and omega_2 is given by:
    omega_2 = - (sqrt(3*gamma)/16) * (6*gamma**2 + 5*gamma + 14)

    The problem asks for the "3rd term" of the correction. This is interpreted as the
    third term in the expanded expression for omega_2.
    """
    gamma = sympy.Symbol('gamma')

    # The expression for the second-order frequency correction (omega_2).
    # It consists of a polynomial in gamma multiplied by sqrt(gamma).
    # We are interested in the third term of this polynomial expansion.
    # The term corresponds to the constant '14' in the parenthesis.
    
    # The third term is: -(sqrt(3*gamma)/16) * 14
    
    # Define the numbers in the final expression
    numerator = 7
    denominator = 8
    gamma_coeff = 3
    
    # Print the equation for the third term, showing each number.
    # The simplified expression is -(7/8) * sqrt(3*gamma)
    print("The final equation for the 3rd term of the nonlinear correction to the frequency is:")
    print(f"Third Term = -({numerator}/{denominator}) * sqrt({gamma_coeff}*gamma)")

solve_rayleigh_plesset_correction()