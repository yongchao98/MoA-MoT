import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def has_prime_factor_4k1(n):
    """Checks if a number has a prime factor of the form 4k+1."""
    # Handle factor 2
    if n % 2 == 0:
        while n % 2 == 0:
            n //= 2
    
    # Handle odd factors
    p = 3
    while p * p <= n:
        if n % p == 0:
            if p % 4 == 1:
                return True
            while n % p == 0:
                n //= p
        p += 2
        
    # Handle the case where n is a prime number itself
    if n > 1 and n % 4 == 1:
        return True
        
    return False

def find_subset():
    """
    Finds the subset of integers t for which the number of t-omino tilings
    of an n x n grid is always even.
    The condition is that t is not a perfect square and has no prime factor of the form 4k+1.
    """
    initial_set = [2, 3, 4, 5, 7, 9, 15]
    result_subset = []
    
    for t in initial_set:
        if not is_perfect_square(t) and not has_prime_factor_4k1(t):
            result_subset.append(t)
            
    print("The subset of integers is:")
    for number in result_subset:
        print(number)

if __name__ == "__main__":
    find_subset()