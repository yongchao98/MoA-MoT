import numpy as np

def solve():
    """
    This script calculates the final answer based on the derived asymptotic formula for l(a).
    """
    n = 1000000000
    
    # The first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    log_n = np.log(n)
    sum_1_over_a = sum(1.0/p for p in primes)
    sum_ln_a_over_a = sum(np.log(p)/p for p in primes)
    
    # Calculate the total sum S
    total_l = -8 * log_n * sum_1_over_a - 16 * sum_ln_a_over_a
    
    # Calculate the floor of the sum
    final_answer = np.floor(total_l)

    print("The final result is the floor of the sum of l(a_i) for the first 10 primes.")
    print("Based on our derivation, the sum is calculated as follows:")
    print(f"n = {n}")
    print(f"a_i = {primes}")
    print(f"ln(n) = {log_n:.5f}")
    print(f"Sum(1/a_i) = {sum_1_over_a:.5f}")
    print(f"Sum(ln(a_i)/a_i) = {sum_ln_a_over_a:.5f}")
    print(f"Sum l(a_i) = -8 * {log_n:.5f} * {sum_1_over_a:.5f} - 16 * {sum_ln_a_over_a:.5f}")
    print(f"Value = {total_l:.5f}")
    print(f"The floor of the value is {int(final_answer)}.")

solve()