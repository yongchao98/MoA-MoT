import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad
import math

def solve():
    """
    This function solves the problem by finding the optimal F, calculating the
    minimized probability of A winning, and then finding the floor of its reciprocal.
    """

    # Step 1: Find the optimal value of F by solving the cubic equation
    # 17*F^3 - 3*F - 2 = 0. We search for a root between 0 and 1.
    def cubic_equation(F):
        return 17 * F**3 - 3 * F - 2

    try:
        sol = root_scalar(cubic_equation, bracket=[0.1, 1.0])
        F_opt = sol.root
    except ValueError:
        print("Could not find the root for F in the given interval.")
        return

    # Step 2: Calculate the minimized probability of A winning for this optimal F.
    # The probability is given by P(F) = F^2/4 + (2/pi) * Integral.
    # The integrand is f(r, F) = r * arccos(sqrt(F(2r-F))/r).
    # The integration is from F/2 to 1.
    def integrand(r, F):
        # The argument of arccos must be between -1 and 1. For r >= F/2,
        # it can be shown that F*(2r-F) <= r^2, so the argument is always valid.
        # We clip it to handle potential floating point inaccuracies.
        arg = np.sqrt(F * (2 * r - F)) / r
        arg = np.clip(arg, -1.0, 1.0)
        return r * np.arccos(arg)

    # Perform the numerical integration
    integral_part, _ = quad(integrand, F_opt / 2, 1, args=(F_opt,))

    # Calculate the total probability
    prob_A_wins = F_opt**2 / 4 + (2 / np.pi) * integral_part

    # Step 3: Compute the final value as requested by the problem.
    final_answer = math.floor(1 / prob_A_wins)
    
    # As requested, output the number from the final equation floor(1/P) = N.
    print(final_answer)

solve()