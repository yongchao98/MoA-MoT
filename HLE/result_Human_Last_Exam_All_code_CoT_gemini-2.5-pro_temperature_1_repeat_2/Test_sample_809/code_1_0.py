import math
from scipy.optimize import fsolve
from scipy.stats import binom

def solve_and_calculate_ev():
    """
    This function follows the step-by-step plan to calculate the expected value.
    """
    # Step 1: Find h
    # The equation for P(f) is P(f) = x + (1-x)P(f)^2.
    # Rearranging gives (1-x)P(f)^2 - P(f) + x = 0.
    # This is a quadratic equation in P(f). For the probability to be guaranteed (P(f)=1),
    # the two roots of the quadratic equation must converge to 1. This happens when the
    # discriminant is zero.
    # Discriminant: b^2 - 4ac = (-1)^2 - 4(1-x)(x) = 1 - 4x + 4x^2 = (1-2x)^2.
    # Setting the discriminant to zero: (1-2x)^2 = 0 => 1-2x = 0 => x = 0.5.
    # So, h is the smallest value of x for which P(f)=1, which is 0.5.
    h = 0.5

    # Step 2: Find P(g)
    # For g, let P1, P2, P3 be the probabilities of the three paths.
    # P1 = h = 0.5
    # P1 + P2 + P3 = 1 => 0.5 + P2 + P3 = 0.5 => P2 + P3 = 0.5
    # P2 = 2 * P3
    # Substituting: 2*P3 + P3 = 0.5 => 3*P3 = 0.5 => P3 = 0.5 / 3 = 1/6
    p3 = 1/6
    # The recursive equation for P(g) is: P(g) = P1*1 + P2*0 + P3*P(g)^6
    # P(g) = 0.5 + (1/6) * P(g)^6
    # Let z = P(g). We need to solve the equation: z^6 - 6z + 3 = 0.
    # We define a function to find its root numerically.
    g_equation = lambda z: z**6 - 6*z + 3
    # The initial guess can be 0.5 as it's the dominant probability term.
    p_g = fsolve(g_equation, 0.5)[0]

    # Step 3: Find P(k)
    # k is a chain of 4 instances of g.
    # P(k) = P(g)^4
    p_k = p_g**4

    # Step 4: Calculate the Expected Value
    # The number of successful traversals S follows a binomial distribution B(n=100, p=p_k).
    n = 100
    p = p_k
    # The opponent wins if S < 6, which means S is 0, 1, 2, 3, 4, or 5.
    # This is the cumulative probability P(S <= 5).
    prob_lose = binom.cdf(k=5, n=n, p=p)
    
    # EV = P(win)*($1) + P(lose)*(-$1) = (1 - prob_lose) - prob_lose = 1 - 2*prob_lose
    expected_value = 1 - 2 * prob_lose

    # Output the results
    print(f"Step 1: The value of h is {h}.")
    print(f"Step 2: The probability of traversing g, P(g), is the root of z^6 - 6z + 3 = 0, which is approximately {p_g:.6f}.")
    print(f"Step 3: The probability of traversing k, P(k), is P(g)^4, which is {p_k:.6f}.")
    print(f"Step 4: The probability of you losing the bet is P(S < 6) = P(S <= 5), which is {prob_lose:.6f}.")
    print("\nThe final equation for the expected value is:")
    print(f"EV = 1 - 2 * P(lose) = 1 - 2 * {prob_lose:.6f} = {expected_value:.6f}")
    print(f"\nYour expected value of playing one game is ${expected_value:.2f}.")

solve_and_calculate_ev()