import math

def solve():
    """
    This function calculates the required value based on the large-n asymptotic
    analysis of the function l(a).
    """
    
    # The dimension n is very large.
    n = 1000000000
    
    # The first 10 prime numbers for a_i.
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    # Based on the derivation, l(a) can be accurately approximated by -n/a.
    # We calculate the sum of l(a_i) over the first 10 primes.
    # The final equation is: floor( l(a_1) + l(a_2) + ... + l(a_10) )
    # where each l(a_i) is approximately -n / a_i.
    
    total_sum = 0
    for p in primes:
        # This is each number in the final equation's sum
        term = -n / p
        total_sum += term
        
    # The final step is to take the floor of the total sum.
    result = math.floor(total_sum)
    
    print(result)

solve()