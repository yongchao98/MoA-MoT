import math

def calculate_probability():
    """
    Calculates the probability that the matrix is diagonalisable.

    The probability is given by the formula:
    P(diag) = 1 - (e^-1 * sum_num / sum_den)
    where:
    sum_num = sum_{k=1 to inf} 1/(k!)^3
    sum_den = sum_{k=0 to inf} 1/(k!)^2
    """
    
    # We approximate the infinite sums by summing up to a large integer N.
    # The factorial terms ensure rapid convergence. N=20 is sufficient for high precision.
    N = 30
    
    sum_numerator = 0.0
    sum_denominator = 0.0
    
    for k in range(N):
        fact_k = math.factorial(k)
        
        # Denominator sum starts from k=0
        sum_denominator += 1 / (fact_k ** 2)
        
        # Numerator sum starts from k=1
        if k > 0:
            sum_numerator += 1 / (fact_k ** 3)
            
    # Calculate e^-1
    e_inv = math.exp(-1)
    
    # Calculate pi, the probability of the matrix being non-diagonalisable
    prob_non_diag = e_inv * sum_numerator / sum_denominator
    
    # The probability of being diagonalisable is 1 - pi
    prob_diag = 1 - prob_non_diag
    
    print("This script calculates the probability that the matrix is diagonalisable.")
    print("The final probability is calculated using the formula:")
    print("P(diagonalisable) = 1 - e^(-1) * (sum_k=1^inf 1/(k!)^3) / (sum_k=0^inf 1/(k!)^2)\n")
    
    print("--- Intermediate Values ---")
    print(f"Value of e^(-1): {e_inv}")
    print(f"Value of the numerator sum (sum_k=1^inf 1/(k!)^3): {sum_numerator}")
    print(f"Value of the denominator sum (sum_k=0^inf 1/(k!)^2): {sum_denominator}")
    print(f"Probability of being non-diagonalisable (pi): {prob_non_diag}\n")
    
    print("--- Final Answer ---")
    print(f"The probability that the matrix is diagonalisable is: {prob_diag}")
    
    # Output the final answer in the required format
    print(f"\n<<<{prob_diag}>>>")

if __name__ == '__main__':
    calculate_probability()