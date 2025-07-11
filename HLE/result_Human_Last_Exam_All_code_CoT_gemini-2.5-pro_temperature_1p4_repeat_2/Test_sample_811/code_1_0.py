import math

def solve_probability():
    """
    Calculates the probability that the matrix is diagonalizable by numerically
    evaluating the derived formula.
    """
    # The sums converge very rapidly due to the factorial term in the Poisson
    # PMF. Truncating at k=30 is more than sufficient for high accuracy.
    MAX_K = 30

    # Pre-calculate Poisson(1) probabilities p_k = e^(-1) / k!
    p = [0.0] * (MAX_K + 1)
    exp_neg_1 = math.exp(-1)
    for k in range(MAX_K + 1):
        p[k] = exp_neg_1 / math.factorial(k)

    # Calculate the numerator sum for the non-diagonalizable probability:
    # N = sum_{k=1 to inf} p_k^3 / (1 + p_k)
    numerator_sum = 0
    for k in range(1, MAX_K + 1):
        # We start from k=1 because the condition for non-diagonalizability is X_N > 0.
        term = p[k]**3 / (1 + p[k])
        numerator_sum += term

    # Calculate the denominator sum:
    # D = sum_{j=0 to inf} p_j^2 / (1 + p_j)
    denominator_sum = 0
    for j in range(MAX_K + 1):
        term = p[j]**2 / (1 + p[j])
        denominator_sum += term

    # Calculate the probability of being non-diagonalizable
    # P(non-diag) = N / D
    prob_non_diag = 0
    if denominator_sum > 1e-20: # Avoid division by zero
        prob_non_diag = numerator_sum / denominator_sum

    # The probability of being diagonalizable is 1 - P(non-diag)
    prob_diag = 1 - prob_non_diag

    # Output the final equation with the computed values
    print("The probability that the matrix is diagonalisable is given by the equation:")
    print(f"P(diagonalisable) = 1 - (Numerator / Denominator)")
    print("\nWhere:")
    print(f"Numerator = sum_{{k=1 to inf}} p_k^3 / (1 + p_k) approx {numerator_sum}")
    print(f"Denominator = sum_{{j=0 to inf}} p_j^2 / (1 + p_j) approx {denominator_sum}")
    
    print("\nFinal Calculation:")
    print(f"1 - ( {numerator_sum} / {denominator_sum} ) = {prob_diag}")

solve_probability()
<<<0.8175850125674066>>>