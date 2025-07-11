import math

def solve_problem():
    """
    This function solves the problem by calculating X_0 and then evaluating the final expression.
    """
    # Step 1: Derivation of X_0.
    # The solvability condition for the boundary-value problem is:
    # integral from 0 to T of (exp(15*t)*x_0^6 + exp(21*t)*y_0^4) dt = 0
    # With T = ln(10), this gives the relation:
    # x_0^6 * (10^15 - 1)/15 + y_0^4 * (10^21 - 1)/21 = 0
    # From this, we can express y_0 in terms of x_0: y_0(x_0) = C * x_0^(3/2),
    # where C = (-7/5 * (10^15 - 1) / (10^21 - 1))^(1/4).

    # The integral condition is given as:
    # integral from 0 to X_0 of y_0(x_0) * x_0^(p-1) dx_0 = beta
    # With p=6, this is integral from 0 to X_0 of (C * x_0^(3/2)) * x_0^5 dx_0 = beta
    # Evaluating the integral gives: C * (2/15) * X_0^(15/2) = beta

    # The value of beta is given as:
    # beta = (1/1000) * (2/15) * C * 10^120
    # (where the term raised to 1/4 is C itself).

    # Equating the two expressions for beta gives:
    # C * (2/15) * X_0^(15/2) = (1/1000) * (2/15) * C * 10^120
    # Cancelling common terms, we get:
    # X_0^(15/2) = 10^(-3) * 10^120 = 10^117

    # Solving for X_0:
    # X_0 = (10^117)^(2/15) = 10^(117 * 2 / 15) = 10^(234 / 15) = 10^15.6
    X0_exponent = 15.6
    X0 = 10**X0_exponent

    # Step 2: Calculate the final expression: 10^30 * X_0^2 - 10^30 * X_0 + 10
    c1 = 10**30
    c2 = 10**30
    c3 = 10

    X0_squared = X0**2
    
    result = c1 * X0_squared - c2 * X0 + c3

    # Step 3: Print the numbers in the final equation as requested.
    print("The value of X_0 is 10**15.6")
    print(f"The value of X_0^2 is 10**{2 * X0_exponent}")
    print("\nThe final equation is:")
    # Using scientific notation for clarity
    print(f"{c1:.1e} * {X0_squared:.5e} - {c2:.1e} * {X0:.5e} + {c3} = {result:.5e}")
    
    # The final answer is the numerical result.
    return result

final_result = solve_problem()
# The final answer is requested in a specific format.
# print(f"\nFinal Answer: {final_result}")
# <<<1.584893192461113e+61>>>