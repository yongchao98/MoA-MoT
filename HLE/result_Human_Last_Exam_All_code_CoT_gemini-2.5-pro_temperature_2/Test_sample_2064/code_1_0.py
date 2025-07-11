import math

def solve():
    """
    This function calculates the final value based on the derived asymptotic formula.
    """
    n = 1000000000
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    log_n = math.log(n)
    
    sum_inv_a = 0
    for p in primes:
        sum_inv_a += 1/p
        
    sum_log_a_div_a = 0
    for p in primes:
        sum_log_a_div_a += math.log(p) / p
        
    # The derived expression for the sum
    total_sum = -log_n * sum_inv_a - 2 * sum_log_a_div_a
    
    # We are asked for the floor of the sum.
    result = math.floor(total_sum)
    
    print("The final result is an integer. The calculation is as follows:")
    print(f"n = {n}")
    print(f"The first 10 prime numbers (a_i) are: {primes}")
    print(f"ln(n) = {log_n}")
    print(f"Sum of 1/a_i = {sum_inv_a}")
    print(f"Sum of ln(a_i)/a_i = {sum_log_a_div_a}")
    print(f"Total sum = -({log_n}) * ({sum_inv_a}) - 2 * ({sum_log_a_div_a})")
    print(f"Total sum = {-log_n * sum_inv_a} + {-2 * sum_log_a_div_a}")
    print(f"Total sum = {total_sum}")
    print(f"The floor of the total sum is {result}")

solve()