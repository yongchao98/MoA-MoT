import math

def solve_probability():
    """
    Calculates the probability that the matrix is diagonalisable.
    
    The probability is given by the formula:
    P(diag) = 1 - e^(-1) * S3 / S2
    where:
    S2 = sum_{k=0 to inf} 1/(k!)^2
    S3 = sum_{k=1 to inf} 1/(k!)^3
    """
    
    s2_sum = 0.0
    s3_sum = 0.0
    
    # The sums converge very quickly. A loop up to 30 is more than sufficient
    # for standard float precision.
    for k in range(30):
        fact_k = math.factorial(k)
        
        # Calculate S2 term
        term_s2 = 1.0 / (fact_k ** 2)
        s2_sum += term_s2
        
        # Calculate S3 term (sum starts from k=1)
        if k >= 1:
            term_s3 = 1.0 / (fact_k ** 3)
            s3_sum += term_s3

    e_inv = math.exp(-1)
    
    # Probability of the matrix not being diagonalisable
    prob_not_diag = e_inv * s3_sum / s2_sum
    
    # The final probability
    prob_diag = 1 - prob_not_diag

    print("The probability is calculated using the formula: P = 1 - e^(-1) * (S3 / S2)")
    print("where the computed values are:")
    print(f"S2 = sum_{{k=0 to inf}} 1/(k!)^2 = {s2_sum}")
    print(f"S3 = sum_{{k=1 to inf}} 1/(k!)^3 = {s3_sum}")
    print(f"e^(-1) = {e_inv}")
    print("\nPlugging in the numbers:")
    print(f"P = 1 - {e_inv} * ({s3_sum} / {s2_sum})")
    print(f"P = 1 - {prob_not_diag}")
    print(f"P = {prob_diag}")

solve_probability()