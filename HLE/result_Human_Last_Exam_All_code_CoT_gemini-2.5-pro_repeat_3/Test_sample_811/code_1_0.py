import math

def calculate_probability():
    """
    Calculates the probability that the given matrix is diagonalisable.
    """
    # The sums converge very quickly. A limit of 30 for k is more than enough.
    LIMIT = 30

    # Let u_k = P(X=k) = e^-1 / k!
    # Let p be the probability that two consecutive variables are equal.
    # p = sum_{k=0 to inf} (u_k)^2 = e^-2 * sum_{k=0 to inf} 1/(k!)^2
    sum_p_terms = 0.0
    for k in range(LIMIT):
        fact_k = math.factorial(k)
        sum_p_terms += 1.0 / (fact_k * fact_k)
    p = math.exp(-2) * sum_p_terms

    # Let q_star be the probability that the first match is a non-zero triple.
    # q_star = sum_{k=1 to inf} (u_k)^3 = e^-3 * sum_{k=1 to inf} 1/(k!)^3
    sum_q_star_terms = 0.0
    for k in range(1, LIMIT):
        fact_k = math.factorial(k)
        sum_q_star_terms += 1.0 / (fact_k * fact_k * fact_k)
    q_star = math.exp(-3) * sum_q_star_terms

    # The probability of being non-diagonalisable is pi = q_star / p
    pi = q_star / p

    # The probability of being diagonalisable is 1 - pi
    prob_diag = 1 - pi

    print("The probability is computed as 1 - P(not diagonalisable).")
    print("P(not diagonalisable) = q* / p, where:")
    print(f"p = sum_{k=0 to inf} [P(X=k)]^2 = {p}")
    print(f"q* = sum_{k=1 to inf} [P(X=k)]^3 = {q_star}")
    print("\nThe final equation with the computed values is:")
    print(f"1 - ( {q_star} / {p} ) = {prob_diag}")

calculate_probability()