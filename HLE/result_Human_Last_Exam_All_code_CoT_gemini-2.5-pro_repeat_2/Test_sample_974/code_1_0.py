def find_prime_divisors(n):
    """
    Finds all unique prime divisors of a given integer n.
    """
    factors = set()
    # Check for divisibility by 2
    if n % 2 == 0:
        factors.add(2)
        while n % 2 == 0:
            n //= 2
    # Check for odd divisors from 3 upwards
    d = 3
    while d * d <= n:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
        d += 2
    # If n is still greater than 1, it must be a prime factor
    if n > 1:
        factors.add(n)
    return sorted(list(factors))

def solve_problem():
    """
    Finds prime divisors p of q such that the number of elements of order p
    in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347
    
    # The condition for the number of elements to be equal is p >= 5,
    # where p is a prime divisor of q.
    
    # Find all prime divisors of q
    prime_divisors_of_q = find_prime_divisors(q)
    
    # Filter the divisors based on the condition p >= 5
    result_primes = []
    for p in prime_divisors_of_q:
        if p >= 5:
            result_primes.append(p)
            
    # The problem asks for an equation. The equation for p is simply that it is one of the found primes.
    print(f"The number q is {q}.")
    print("The number of elements of order p in PSL(3, q^2) and PSL(4, q) are equal if p is a prime divisor of q and p >= 5.")
    print("The prime divisors of q satisfying this condition are:")
    if not result_primes:
        print("None")
    else:
        for p in result_primes:
            # The prompt asks to "output each number in the final equation".
            # We interpret this as showing the resulting values for p.
            print(f"p = {p}")

solve_problem()
