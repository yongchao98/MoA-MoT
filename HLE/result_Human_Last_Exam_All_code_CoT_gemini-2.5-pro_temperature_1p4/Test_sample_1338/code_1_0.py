import math

def is_prime_power(n):
    """
    Checks if a number n is a prime power (of the form p^k for prime p, k>=1).
    """
    if n <= 1:
        return False
    
    # Find the prime factorization of n.
    # A number is a prime power if it has exactly one distinct prime factor.
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
        
    return len(factors) == 1

def find_subset():
    """
    Finds the subset of t values based on the prime power hypothesis.
    """
    t_values = [2, 3, 4, 5, 7, 9, 15]
    result_set = []

    for t in t_values:
        if not is_prime_power(t):
            result_set.append(t)
    
    print("The subset of integers for which the statement is true is:")
    # The problem asks to output each number in the final equation.
    # We will print each number that forms our resulting set.
    for number in result_set:
        print(number)

if __name__ == '__main__':
    find_subset()
