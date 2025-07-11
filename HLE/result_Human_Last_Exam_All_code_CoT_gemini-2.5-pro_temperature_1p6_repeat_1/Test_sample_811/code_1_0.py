import math

def calculate_probability():
    """
    Calculates the probability that the matrix is diagonalisable.
    
    The probability is given by the formula:
    P(diag) = 1 - (sum_{k=1 to inf} p_k^3) / (sum_{k=0 to inf} p_k^2)
    where p_k = exp(-1) / k! is the Poisson(1) probability mass function.
    
    We approximate the infinite sums by summing up to a large number of terms.
    """
    
    # Use a sufficient number of terms for the sums to converge
    max_k = 30
    
    # Calculate p_k = exp(-1) / k! for k from 0 to max_k
    p_values = [math.exp(-1) / math.factorial(k) for k in range(max_k + 1)]
    
    # Calculate the denominator: S2 = sum_{k=0 to inf} p_k^2
    sum_pk_2 = sum(p**2 for p in p_values)
    
    # Calculate the numerator: S3_prime = sum_{k=1 to inf} p_k^3
    # We sum from k=1 onwards, so we skip the first p_value (for k=0)
    sum_pk_3_from_1 = sum(p**3 for p in p_values[1:])
    
    # Calculate the probability of being non-diagonalisable
    prob_not_diag = sum_pk_3_from_1 / sum_pk_2
    
    # The probability of being diagonalisable is 1 - prob_not_diag
    prob_diag = 1 - prob_not_diag
    
    print("The final equation for the probability of being diagonalisable is:")
    print("P(Diagonalisable) = 1 - (Numerator / Denominator)")
    print(f"Denominator = sum_{k=0 to inf}(p_k^2) = {sum_pk_2}")
    print(f"Numerator = sum_{k=1 to inf}(p_k^3) = {sum_pk_3_from_1}")
    print(f"P(Diagonalisable) = 1 - ({sum_pk_3_from_1} / {sum_pk_2}) = {prob_diag}")

calculate_probability()