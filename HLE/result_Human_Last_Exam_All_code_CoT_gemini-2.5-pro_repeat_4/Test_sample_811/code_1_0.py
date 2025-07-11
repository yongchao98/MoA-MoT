import math

def solve_probability():
    """
    Calculates the probability that the given matrix is diagonalisable.
    
    The probability P is given by the formula:
    P = 1 - (e^-1 * (S3 - 1)) / S2
    where:
    S2 = sum_{k=0 to inf} 1 / (k!)^2
    S3 = sum_{k=0 to inf} 1 / (k!)^3
    """
    
    # Plan:
    # 1. Calculate S2 and S3 by summing a sufficient number of terms.
    #    The series converge very rapidly, so 30 terms will provide high precision.
    # 2. Calculate e^-1.
    # 3. Substitute these values into the formula to get the final probability.
    # 4. Print the intermediate values and the final result as requested.
    
    n_terms = 30
    s2 = 0.0
    s3 = 0.0
    
    # Iteratively compute sums for S2 and S3
    for k in range(n_terms):
        try:
            # math.factorial is efficient and sufficient for this range
            fact_k = math.factorial(k)
            # The term 1/(k!) is calculated once per iteration
            inv_fact_k = 1 / fact_k
            
            s2 += inv_fact_k ** 2
            s3 += inv_fact_k ** 3
        except OverflowError:
            # For very large k, factorial will overflow, but the terms
            # are already negligibly small, so we can stop.
            break
            
    e_inv = math.exp(-1)
    
    # The derived formula for the probability of the matrix being diagonalisable
    final_prob = 1 - e_inv * (s3 - 1) / s2
    
    print("This script calculates the probability that the specified matrix is diagonalisable.")
    print("The final probability is calculated using the formula:")
    print("P = 1 - (exp(-1) * (S3 - 1)) / S2")
    print("where S2 = sum_k(1/(k!)^2) and S3 = sum_k(1/(k!)^3).")
    print("\nCalculating the numerical values...")
    print("-" * 50)
    print(f"Value of S2 (sum of 1/(k!)^2): {s2}")
    print(f"Value of S3 (sum of 1/(k!)^3): {s3}")
    print(f"Value of exp(-1): {e_inv}")
    print("-" * 50)
    print("The final equation with the computed values is:")
    print(f"P = 1 - ({e_inv} * ({s3} - 1)) / {s2}")
    print("\nFinal Result:")
    print(f"The probability is: {final_prob}")

solve_probability()