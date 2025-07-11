import math

def calculate_probability():
    """
    Calculates the probability that the matrix is diagonalisable.
    
    The probability is given by the formula:
    P(diagonalisable) = 1 - P(not diagonalisable)
    
    where
    
    P(not diagonalisable) = Numerator / Denominator
    
    Numerator = sum_{k=1 to inf} [ p_k^3 / (1 + p_k) ]
    Denominator = sum_{j=0 to inf} [ p_j^2 / (1 + p_j) ]
    
    and p_k is the Poisson(1) probability P(X=k) = exp(-1) / k!
    """
    
    # We can truncate the sum at a sufficiently large number, as the terms decrease very rapidly.
    max_k = 20
    
    # Pre-compute Poisson probabilities p_k = exp(-1)/k!
    try:
        p_values = [math.exp(-1) / math.factorial(k) for k in range(max_k)]
    except (ValueError, OverflowError):
        # Fallback for platforms with smaller limits for factorial
        p_values = []
        for k in range(max_k):
            try:
                p_values.append(math.exp(-1) / math.factorial(k))
            except (ValueError, OverflowError):
                break
    
    # Calculate the numerator: sum from k=1
    numerator_sum = sum(p**3 / (1 + p) for p in p_values[1:])
    
    # Calculate the denominator: sum from j=0
    denominator_sum = sum(p**2 / (1 + p) for p in p_values)
    
    # Calculate the probability of not being diagonalisable
    prob_not_diag = numerator_sum / denominator_sum
    
    # The probability of being diagonalisable is 1 minus that.
    prob_diag = 1 - prob_not_diag

    print("The final probability is derived from the expression: 1 - Numerator / Denominator")
    print(f"Numerator = {numerator_sum}")
    print(f"Denominator = {denominator_sum}")
    print(f"Probability of not being diagonalisable = {prob_not_diag}")
    print(f"The probability that the matrix is diagonalisable is: {prob_diag}")

calculate_probability()