import math

def find_prime_factors(n):
    """
    This function finds the unique prime factors of a given integer n.
    It returns a sorted list of the prime factors.
    """
    factors = set()
    # Handle the factor 2 separately
    while n % 2 == 0:
        factors.add(2)
        n //= 2
    
    # Iterate for odd factors from 3 up to the square root of n
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.add(i)
            n //= i
        i += 2
        
    # If n is a prime number greater than 2, it remains after the loop
    if n > 1:
        factors.add(n)
        
    return sorted(list(factors))

def solve_group_theory_problem():
    """
    Solves the main problem by finding the required prime divisors.
    
    The problem reduces to finding all prime divisors p of q such that p >= 5.
    As q is not divisible by 2, 3, or 5, all its prime factors will meet this condition.
    Therefore, the task is simply to list all prime factors of q.
    """
    q = 12740347
    
    prime_divisors = find_prime_factors(q)
    
    print(f"The prime divisors p of q={q} for which the number of elements of order p in both PSL(3,q^2) and PSL(4,q) are equal are:")
    for p in prime_divisors:
        # The condition is that p >= 5, which holds for all prime factors of this specific q.
        print(p)

solve_group_theory_problem()