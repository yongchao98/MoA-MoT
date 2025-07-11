import math
from decimal import Decimal, getcontext

def solve_problem():
    """
    This function solves the problem by following these analytical steps:

    1.  The solvability condition for the given boundary-value problem establishes a relationship between the initial values x_0 and y_0:
        x_0**p * Integral_1 + y_0**q * Integral_2 = -alpha
        Given p=6, q=4, A=3, B=2A=6, T=ln(10), alpha=0, the condition becomes:
        x_0**6 * (10**15 - 1)/15 + y_0**4 * (10**21 - 1)/21 = 0
        This yields y_0**4 = C * x_0**6, where C is a constant.
        Thus, y_0(x_0) = C**(1/4) * x_0**(6/4) = C**(1/4) * x_0**(3/2).

    2.  The problem provides an integral equation for a point X_0:
        integral from 0 to X_0 of y_0(x_0) * x_0**(p-1) dx_0 = beta
        Substituting y_0(x_0) and p=6:
        integral from 0 to X_0 of (C**(1/4) * x_0**(3/2)) * x_0**5 dx_0 = beta
        C**(1/4) * integral from 0 to X_0 of x_0**(13/2) dx_0 = beta
        C**(1/4) * [x_0**(15/2) / (15/2)] evaluated at X_0 = beta
        This gives: C**(1/4) * (2/15) * X_0**(15/2) = beta

    3.  The value of beta is given as:
        beta = (1/1000) * (2/15) * C**(1/4) * 10**120
        where the term under the 1/4 power is indeed the constant C.

    4.  Equating the two expressions for beta, the terms C**(1/4) and 2/15 cancel out:
        X_0**(15/2) = (1/1000) * 10**120
        X_0**(15/2) = 10**(-3) * 10**120 = 10**117

    5.  Solving for X_0:
        X_0 = (10**117)**(2/15) = 10**(117 * 2 / 15) = 10**(15.6)
    """

    # Set precision for Decimal calculations to handle large numbers.
    getcontext().prec = 100

    # From the derivation, the exponent for X_0 is 15.6
    exp_X0 = Decimal('15.6')
    X0 = Decimal(10) ** exp_X0

    # The final expression to calculate is 10^30 * X_0^2 - 10^30 * X_0 + 10
    coeff = Decimal(10)**Decimal(30)
    
    # Calculate each term of the final equation
    X0_sq = X0 ** 2
    term1 = coeff * X0_sq
    term2 = coeff * X0
    term3 = Decimal(10)

    # Calculate the final result
    result = term1 - term2 + term3

    print("The final equation is: Term1 - Term2 + Term3")
    print(f"where X_0 = 10**{exp_X0}")
    print(f"Term1 = 10**30 * X_0**2 = 10**30 * (10**{exp_X0})**2 = 10**(30 + 2*{exp_X0}) = 10**{30 + 2*exp_X0}")
    print(f"Term2 = 10**30 * X_0 = 10**30 * 10**{exp_X0} = 10**(30 + {exp_X0}) = 10**{30 + exp_X0}")
    print(f"Term3 = {term3}")
    print("\nCalculating the numerical value:")
    # The result is a very large number, which we print in scientific notation.
    print(f"Result = {term1:.4e} - {term2:.4e} + {term3}")
    print(f"Final Answer = {result:.4e}")

solve_problem()