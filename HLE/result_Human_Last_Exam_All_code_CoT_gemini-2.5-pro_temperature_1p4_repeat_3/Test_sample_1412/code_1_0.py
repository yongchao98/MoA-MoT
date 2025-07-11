import math

def get_divisors(n):
    """Returns all divisors of a number n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def euler_totient(n):
    """Calculates Euler's totient function phi(n)."""
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
    Calculates the number of non-isomorphic, connected, 3-regular,
    adjustable graphs with 2000 vertices that have a perfect matching.
    """
    n = 1000
    k = 2  # Number of choices (colors) for each edge

    # Calculate sum over rotations for Polya Enumeration Theorem
    # Sum_{d|n} phi(d) * k^(n/d)
    divs = get_divisors(n)
    rotation_sum = 0
    for d in divs:
        phi_d = euler_totient(d)
        power_term = pow(k, n // d)
        rotation_sum += phi_d * power_term

    # Handle reflections
    # For n even, there are n/2 reflections through vertices and n/2 through edges.
    # Term is (n/2)*k^(n/2) + (n/2)*k^(n/2+1)
    if n % 2 == 0:
        reflection_sum = (n // 2) * pow(k, n // 2) + (n // 2) * pow(k, (n // 2) + 1)
    else: # n is odd
        reflection_sum = n * pow(k, (n + 1) // 2)

    # Total number of non-isomorphic graphs is the sum divided by group size (2n)
    total_sum = rotation_sum + reflection_sum
    num_graphs = total_sum // (2 * n)
    
    print("The number of non-isomorphic graphs is:")
    print(num_graphs)

solve()