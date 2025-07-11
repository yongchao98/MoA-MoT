import math

def solve():
    """
    This function calculates the cardinality of Theta^{-1}(lambda) for m=3.
    """
    m = 3
    
    # Calculate n for m=3
    # n = sum_{k=1 to m} k*(m+1-k)
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))
    
    # For m=3, the partition lambda of n=10 is (3^1, 2^2, 1^3)
    # This means one part of size 3, two parts of size 2, and three parts of size 1.
    partition_counts = {
        3: 1,
        2: 2,
        1: 3,
    }

    # Calculate the size of the centralizer, z_lambda
    # z_lambda = product_{k} k^{i_k} * i_k!
    z_lambda = 1
    for k, i_k in partition_counts.items():
        z_lambda *= (k**i_k) * math.factorial(i_k)
        
    # Calculate n!
    n_factorial = math.factorial(n)
    
    # The cardinality is given by the formula (n!)^2 / z_lambda
    # We use integer arithmetic to maintain precision
    numerator = n_factorial * n_factorial
    cardinality = numerator // z_lambda
    
    # Output the result as an equation
    print(f"For m = {m}, we have n = {n}.")
    print(f"The partition lambda is (3^1, 2^2, 1^3).")
    print(f"The size of the centralizer z_lambda is {z_lambda}.")
    print(f"The cardinality is (n!)^2 / z_lambda = ({n_factorial})^2 / {z_lambda} = {cardinality}.")
    
    # Get the first 40 digits of the result
    result_str = str(cardinality)
    first_40_digits = result_str[:40]
    
    print("\nThe first 40 digits of the cardinality are:")
    print(first_40_digits)

solve()