import math

def get_prime_divisors(n):
    """
    Returns a list of unique prime divisors of an integer n.
    """
    factors = set()
    # Handle the factor 2
    if n % 2 == 0:
        factors.add(2)
        while n % 2 == 0:
            n //= 2
            
    # Handle odd factors
    d = 3
    # We only need to check up to the square root of n
    limit = int(math.sqrt(n))
    while d <= limit:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
            # Update the limit, as n has been reduced
            limit = int(math.sqrt(n))
        d += 2
        
    # If n is still greater than 1 after the loop, n itself is a prime factor
    if n > 1:
        factors.add(n)
        
    return sorted(list(factors))

def solve_problem():
    """
    Solves the problem by finding the prime divisors for which the condition holds.
    """
    q = 12740347

    # As derived in the explanation, the number of elements of order p in both
    # PSL(3, q^2) and PSL(4, q) is q^12 - 1.
    # The equation is q^12 - 1 = q^12 - 1, which is always true.
    # Therefore, we just need to find all prime divisors of q.
    
    prime_divisors_of_q = get_prime_divisors(q)
    
    print(f"Let q = {q}.")
    print("The number of elements of order p in PSL(3, q^2) is q^12 - 1.")
    print("The number of elements of order p in PSL(4, q) is q^12 - 1.")
    print("The final equation is q^12 - 1 = q^12 - 1, which holds for any q > 1.")
    print("This condition is satisfied for all prime divisors p of q.")
    print(f"\nThe prime divisors of q = {q} are:")
    print(prime_divisors_of_q)

solve_problem()
