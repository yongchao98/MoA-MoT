import math

def calculate_probability():
    """
    Calculates the probability that the given matrix is diagonalizable.
    """
    # Set a truncation limit for the infinite sums.
    # The terms decrease very quickly, so K=25 is more than sufficient for high precision.
    K = 25
    
    # Pre-compute Poisson probabilities p_k = P(X=k) for X ~ Poisson(1)
    p = [0.0] * K
    for k in range(K):
        try:
            p[k] = math.exp(-1) / math.factorial(k)
        except (ValueError, OverflowError):
            # factorial(k) can become very large
            p[k] = 0

    # Calculate the numerator sum for P(not diag)
    # This is Sum_{k=1 to inf} [ p_k^3 / (1 + p_k) ]
    numerator_sum = 0
    for k in range(1, K):
        if p[k] > 0:
            numerator_sum += (p[k]**3) / (1 + p[k])
            
    # Calculate the denominator sum for P(not diag)
    # This is Sum_{j=0 to inf} [ p_j^2 / (1 + p_j) ]
    denominator_sum = 0
    for j in range(K):
        if p[j] > 0:
            denominator_sum += (p[j]**2) / (1 + p[j])

    # Calculate the probability of the matrix not being diagonalizable
    if denominator_sum == 0:
        prob_not_diag = 0
    else:
        prob_not_diag = numerator_sum / denominator_sum

    # The probability of being diagonalizable is 1 minus the above
    prob_diag = 1 - prob_not_diag

    print(f"The numerator of the P(not diag) expression is: {numerator_sum}")
    print(f"The denominator of the P(not diag) expression is: {denominator_sum}")
    print(f"The probability that the matrix is NOT diagonalizable is: {prob_not_diag}")
    print(f"The probability that the matrix IS diagonalizable is: {prob_diag}")

calculate_probability()