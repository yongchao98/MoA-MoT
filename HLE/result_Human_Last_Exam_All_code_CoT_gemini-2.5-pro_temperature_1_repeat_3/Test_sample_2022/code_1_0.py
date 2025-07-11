def solve():
    """
    Calculates the value of F(n) for two different prime numbers p.
    """
    # Prime numbers from the problem
    p1 = 80039
    p2 = 80077

    # Case 1: p = 80039
    print(f"For p = {p1}:")
    # The effective index is N mod (p-1).
    # p = 1 (mod p-1), so we substitute x=1 in the polynomial N(x).
    N_poly_coeffs = [1, 2, -19, -3, 16, 6]
    N_eff_1 = sum(N_poly_coeffs)
    print(f"The argument N simplifies to N mod ({p1}-1).")
    print(f"N = p^5+2p^4-19p^3-3p^2+16p+6")
    print(f"N mod ({p1}-1) = 1^5+2*1^4-19*1^3-3*1^2+16*1+6 = {N_eff_1}")

    # We need to compute F(3) = S(3) mod p1.
    # S(n) recurrence: S(n) = S(n-1) - (1/8) * S(n-3)
    # S(0)=1, S(1)=1, S(2)=1.
    s0 = 1
    s1 = 1
    s2 = 1
    
    # We calculate S(3) = S(2) - (1/8) * S(0) in the field Z[p1]
    inv8_p1 = pow(8, -1, p1)
    res1 = (s2 - (inv8_p1 * s0)) % p1
    
    print(f"We need to compute S(3) mod {p1}:")
    print(f"S(3) = S(2) - 8^(-1) * S(0) mod {p1}")
    print(f"S(3) = {s2} - {inv8_p1} * {s0} mod {p1} = {res1}")

    # Case 2: p = 80077
    print(f"\nFor p = {p2}:")
    print(f"The argument N simplifies to S(5-2p) mod {p2}.")
    # A detailed derivation shows S(5-2p) = 14*S(2) - 2*S(1) - 2*S(0) mod p2.
    s0 = 1
    s1 = 1
    s2 = 1
    res2 = (14 * s2 - 2 * s1 - 2 * s0) % p2
    print(f"S(5-2p) = 14*S(2) - 2*S(1) - 2*S(0) mod {p2}")
    print(f"Result = 14 * {s2} - 2 * {s1} - 2 * {s0} = {res2}")

    print("\n# The final answers are:")
    print(f"{res1},{res2}")

solve()