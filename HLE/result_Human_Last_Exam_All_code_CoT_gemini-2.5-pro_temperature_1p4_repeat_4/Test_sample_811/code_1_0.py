import math

def solve_probability():
    """
    Calculates the probability that the given matrix is diagonalizable.

    The probability is given by the formula P = 1 - T/S, where:
    - p_k is the probability mass function for a Poisson(1) distribution.
    - S is the sum of p_k^2 over all k >= 0.
    - T is the sum of p_k^3 over all k >= 1.
    """
    n_terms = 30  # Series converge very fast, 30 terms is sufficient for high precision.
    
    # Pre-calculate exp(-1)
    exp_neg_1 = math.exp(-1)
    
    sum_pk_2 = 0.0
    sum_pk_3 = 0.0
    
    # We can pre-calculate factorials for efficiency
    factorials = [1] * (n_terms)
    for i in range(1, n_terms):
        factorials[i] = factorials[i-1] * i

    # Calculate S and T by summing the series
    for k in range(n_terms):
        p_k = exp_neg_1 / factorials[k]
        
        pk_2 = p_k ** 2
        sum_pk_2 += pk_2
        
        if k >= 1:
            pk_3 = p_k ** 3
            sum_pk_3 += pk_3
            
    S = sum_pk_2
    T = sum_pk_3
    
    probability = 1 - (T / S)
    
    print("The final equation for the probability is 1 - T/S, where T is summed from k=1 and S is summed from k=0.")
    print(f"S = sum(p_k^2) = {S}")
    print(f"T = sum(p_k^3) for k>0 = {T}")
    print(f"Probability = 1 - {T}/{S} = {probability}")

solve_probability()