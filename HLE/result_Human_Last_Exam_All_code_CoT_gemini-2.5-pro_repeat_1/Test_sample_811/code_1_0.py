import math

def solve_probability():
    """
    Calculates the probability that the given matrix is diagonalizable.
    
    The probability is given by the formula:
    P(diag) = 1 - P(not diag)
    where P(not diag) = (sum_{k=1 to inf} p_k^3 / (1 + p_k)) / (sum_{j=0 to inf} p_j^2 / (1 + p_j))
    and p_k is the Poisson(1) probability P(X=k).
    """
    
    # We truncate the infinite sums at a large value K for numerical approximation.
    # p_k becomes extremely small for large k, so K=30 is more than sufficient.
    K = 30
    
    numerator_sum = 0.0
    denominator_sum = 0.0
    
    # Pre-calculate Poisson probabilities to avoid re-computation
    poisson_probs = []
    for i in range(K + 1):
        try:
            # P(X=i) for X ~ Poisson(1) is e^-1 / i!
            prob = math.exp(-1) / math.factorial(i)
            poisson_probs.append(prob)
        except (ValueError, OverflowError):
            # Factorial can get very large, break if we hit limits
            break

    # Calculate the denominator sum: sum_{j=0 to inf} p_j^2 / (1 + p_j)
    for j in range(len(poisson_probs)):
        pj = poisson_probs[j]
        denominator_sum += (pj**2) / (1 + pj)
        
    # Calculate the numerator sum: sum_{k=1 to inf} p_k^3 / (1 + p_k)
    # The sum starts from k=1
    for k in range(1, len(poisson_probs)):
        pk = poisson_probs[k]
        numerator_sum += (pk**3) / (1 + pk)
        
    # Calculate the probability of not being diagonalizable
    prob_not_diag = numerator_sum / denominator_sum
    
    # The probability of being diagonalizable is 1 minus the above
    prob_diag = 1 - prob_not_diag
    
    print("The probability of the matrix being diagonalizable is P(diag) = 1 - P(not diag).")
    print("\nFirst, we approximate the components of P(not diag):")
    print(f"Numerator sum (approximated) = {numerator_sum}")
    print(f"Denominator sum (approximated) = {denominator_sum}")
    
    print("\nThen, we compute P(not diag):")
    print(f"P(not diag) = {numerator_sum} / {denominator_sum} = {prob_not_diag}")
    
    print("\nFinally, we compute the desired probability, P(diag):")
    print(f"1 - {prob_not_diag} = {prob_diag}")

solve_probability()