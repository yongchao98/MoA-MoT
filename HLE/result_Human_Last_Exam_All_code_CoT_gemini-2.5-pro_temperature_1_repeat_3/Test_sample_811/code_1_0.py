import math

def calculate_probability():
    """
    Calculates the probability that the given matrix is diagonalisable.
    """
    # The probability mass function for a Poisson(1) random variable is p_k = e^(-1) / k!
    # We need to compute two sums:
    # 1. Numerator for P(not diag): sum_{k=1 to inf} (p_k)^3
    # 2. Denominator for P(not diag): sum_{j=0 to inf} (p_j)^2
    
    # The terms in the sums decrease very rapidly, so we can approximate the infinite
    # sums by summing a finite number of terms (e.g., up to 30).
    
    # Pre-calculate factorials to speed up computation
    factorials = [1] * 31
    for i in range(1, 31):
        factorials[i] = factorials[i-1] * i
        
    e_inv = 1 / math.e
    
    # Calculate numerator: sum_{k=1 to inf} p_k^3
    # Note that the sum starts from k=1 for the non-diagonalisable case.
    numerator_sum = 0
    for k in range(1, 31):
        p_k = e_inv / factorials[k]
        numerator_sum += p_k**3
        
    # Calculate denominator: sum_{j=0 to inf} p_j^2
    denominator_sum = 0
    for j in range(31):
        p_j = e_inv / factorials[j]
        denominator_sum += p_j**2
        
    # Probability of being non-diagonalisable
    prob_not_diag = numerator_sum / denominator_sum
    
    # The probability of being diagonalisable is 1 minus the probability of not being diagonalisable
    prob_diag = 1 - prob_not_diag
    
    print("The probability of the matrix being diagonalisable is given by the equation:")
    print("P(diagonalisable) = 1 - (sum_{k=1 to inf} p_k^3) / (sum_{j=0 to inf} p_j^2)")
    print("\nCalculating the components:")
    print(f"Numerator sum (sum_{k=1 to inf} p_k^3) = {numerator_sum}")
    print(f"Denominator sum (sum_{j=0 to inf} p_j^2) = {denominator_sum}")
    print(f"P(not diagonalisable) = {numerator_sum} / {denominator_sum} = {prob_not_diag}")
    print(f"P(diagonalisable) = 1 - {prob_not_diag} = {prob_diag}")

calculate_probability()