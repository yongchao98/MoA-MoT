import numpy as np

def solve_ode_approximation():
    """
    This function derives the approximate analytical solution for the given ODE in the large x regime.
    """

    # As derived in the plan, the dominant balance y''' ≈ (y')^4 near the singularity
    # leads to the condition A^3 = -30 for the coefficient in the solution form y(x) = C + A(x0-x)^(1/3).
    A_cubed = -30
    A = np.cbrt(A_cubed)

    # From the ratio of initial conditions y''(0)/y'(0) = 2/3, we find x0 = 1.
    x0 = 1.0

    # From the initial condition y(0)=0, we find C = -A.
    C = -A

    # Round the constant to two decimal places for the final expression.
    coeff = round(C, 2)

    # The power in the expression is n = 1/3.
    power_numerator = 1
    power_denominator = 3

    # Print the final expression with each number explicitly stated.
    print("The analytical expression that approximates the solution for large x is:")
    # The final expression is y(x) = C - A * (x0 - x)^(1/3)
    # y(x) = C * (1 - (x0 - x)^(1/3)), substituting A=-C
    print(f"y(x) = {coeff:.2f} * (1 - ({int(x0)} - x)**({power_numerator}/{power_denominator}))")

solve_ode_approximation()
