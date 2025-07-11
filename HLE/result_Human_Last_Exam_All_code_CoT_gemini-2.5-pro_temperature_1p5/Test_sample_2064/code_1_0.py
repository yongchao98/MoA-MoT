import numpy as np

def solve():
    """
    Calculates the final value based on the derived approximation.
    """
    n = 1_000_000_000
    
    # The first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    # The sum is dominated by the term -(n+1)/a
    # sum_{i=1 to 10} l(a_i) approx -(n+1) * sum_{i=1 to 10} 1/a_i
    
    sum_of_reciprocal_primes = sum(1/p for p in primes)
    
    total_sum = -(n + 1) * sum_of_reciprocal_primes
    
    # Calculate the floor of the result
    final_answer = np.floor(total_sum)
    
    # The final expression is -(n+1) * (1/a_1 + 1/a_2 + ... + 1/a_10)
    print("n = 1000000000")
    # For transparency, let's print each term of the sum of reciprocals
    sum_expr_parts = []
    for p in primes:
        sum_expr_parts.append(f"1/{p}")
    print(f"Sum of reciprocals = {' + '.join(sum_expr_parts)}")
    print(f"Sum of reciprocals value = {sum_of_reciprocal_primes}")

    final_expr = f"floor( -({n}+1) * ({sum_of_reciprocal_primes}) )"
    print(f"Calculation: {final_expr}")
    print(f"Result: {final_answer}")
    
solve()