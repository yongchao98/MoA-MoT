import math

def calculate_probability():
    """
    Calculates the probability that the matrix is diagonalisable.
    The formula is P = 1 - exp(-1) * (sum_{k=1 to inf} 1/(k!)^3) / (sum_{k=0 to inf} 1/(k!)^2).
    We approximate the infinite sums by summing up to a finite limit max_k.
    """
    # The upper limit for the summation. Terms 1/(k!)^n decrease very rapidly.
    max_k = 30

    # sum_denom_terms approximates sum_{k=0 to inf} 1/(k!)^2
    sum_denom_terms = 0.0
    # sum_num_terms approximates sum_{k=1 to inf} 1/(k!)^3
    sum_num_terms = 0.0

    # Iteratively calculate k! to avoid recomputation and large number issues.
    fact_k = 1.0

    # Handle the k=0 term for the denominator sum
    sum_denom_terms += 1.0 / (fact_k * fact_k)

    # Loop from k=1 up to max_k
    for k in range(1, max_k + 1):
        fact_k *= k
        # Stop if factorial becomes too large for float representation
        if fact_k == float('inf'):
            break
        
        # Denominator term: 1/(k!)^2
        term_denom = 1.0 / (fact_k * fact_k)
        sum_denom_terms += term_denom
        
        # Numerator term: 1/(k!)^3
        term_num = term_denom / fact_k
        sum_num_terms += term_num

    # Constant e^-1
    exp_minus_1 = 1 / math.e
    
    # Calculate the probability of the non-diagonalisable event
    prob_non_diag = exp_minus_1 * sum_num_terms / sum_denom_terms
    
    # The probability of being diagonalisable is 1 minus the complement's probability
    prob_diag = 1 - prob_non_diag

    # Output the numbers used in the final calculation, as requested
    print("The final probability is calculated as: 1 - P(not diagonalisable)")
    print(f"P(not diagonalisable) = exp(-1) * (numerator_sum) / (denominator_sum)")
    print(f"P(not diagonalisable) = {exp_minus_1} * {sum_num_terms} / {sum_denom_terms}")
    print(f"P(not diagonalisable) = {prob_non_diag}")
    print("\nResult:")
    print(f"The probability that the matrix is diagonalisable is: {prob_diag}")

calculate_probability()