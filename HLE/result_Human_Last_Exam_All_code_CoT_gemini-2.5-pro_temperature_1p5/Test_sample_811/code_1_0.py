import math

def calculate_probability():
    """
    Calculates the probability that the matrix is diagonalisable.
    """
    # We compute the sums up to a certain number of terms. The series converge very rapidly.
    max_terms = 20
    
    # Numerator part: sum_{k=1 to infinity} 1/(k!)^3
    sum_num = 0.0
    # Denominator part: sum_{j=0 to infinity} 1/(j!)^2
    sum_den = 0.0
    
    for k in range(max_terms):
        fact_k = math.factorial(k)
        
        # Add term to the denominator sum
        term_den = 1.0 / (fact_k ** 2)
        sum_den += term_den
        
        # Add term to the numerator sum (starting from k=1)
        if k >= 1:
            term_num = 1.0 / (fact_k ** 3)
            sum_num += term_num
            
    # Calculate e^-1
    e_inv = math.exp(-1)
    
    # Probability of NOT being diagonalisable
    prob_not_diag = e_inv * sum_num / sum_den
    
    # Probability of being diagonalisable
    prob_diag = 1 - prob_not_diag
    
    print("The formula for the probability is P = 1 - exp(-1) * A / B")
    print("where:")
    print(f"  A = sum_{{k=1 to inf}} 1/(k!)^3 approx {sum_num}")
    print(f"  B = sum_{{k=0 to inf}} 1/(k!)^2 approx {sum_den}")
    print(f"  exp(-1) approx {e_inv}")
    print("\nPlugging in the numbers:")
    print(f"P = 1 - {e_inv} * {sum_num} / {sum_den}")
    print(f"P = 1 - {prob_not_diag}")
    print(f"\nThe final probability is: {prob_diag}")

calculate_probability()