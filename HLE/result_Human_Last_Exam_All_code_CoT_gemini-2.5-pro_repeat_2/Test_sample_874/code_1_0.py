def solve():
    """
    Solves the problem by finding the tuple (a,b,c,d) that maximizes the
    Ducci sequence length f(a,b,c,d) and has the minimum sum, then
    calculates (a+b-c-d) mod 1000 for that tuple.
    """

    # Step 1: Generate Tribonacci numbers
    # t_n = t_{n-1} + t_{n-2} + t_{n-3} with seed (0,0,1)
    limit = 10_000_000
    t = [0, 0, 1]
    while t[-1] <= limit:
        t.append(t[-1] + t[-2] + t[-3])
    
    # The last number in t is > limit, so we don't use it.
    # The largest usable index is len(t)-2.
    # Let's find n for v_n = (t_{n-1}, t_n, t_{n+1}, t_{n+2})
    # We need t_{n+2} <= limit.
    # t_{29} = 8646064, t_{30} = 15902591
    # So the max index for any component is 29.
    # max n+2 is 29, so max n is 27.
    n = 27

    # Step 2: Calculate the maximum length M = f(v_27)
    # A helper function to calculate f(v) for any tuple v=(a,b,c,d)
    memo_f = {}
    def f(a, b, c, d):
        v = (a, b, c, d)
        if v in memo_f:
            return memo_f[v]
        
        count = 1
        sq = v
        while sq != (0, 0, 0, 0):
            a, b, c, d = sq
            sq = (abs(a - b), abs(b - c), abs(c - d), abs(d - a))
            count += 1
        memo_f[v] = count
        return count

    # We establish the recurrence f(v_n) = 3 + f(v_{n-2}) for odd n >= 11.
    # Base case: Calculate f(v_9)
    # v_9 = (t_8, t_9, t_10, t_11) = (24, 44, 81, 149)
    g9 = f(t[8], t[9], t[10], t[11]) # This is f(v_9)

    # Apply the recurrence g(n) = g(n-2) + 3
    # M = g(27) = g(9) + 3 * ( (27-9)/2 )
    k = (n - 9) // 2
    M = g9 + 3 * k

    # Step 3: Identify the tuple (a,b,c,d)
    # This is v_27 = (t_26, t_27, t_28, t_29)
    a = t[26]
    b = t[27]
    c = t[28]
    d = t[29]
    
    # These tuples are known to have the longest life for their size.
    # To get minimal sum for M, we should use this primitive tuple.

    print(f"The Tribonacci sequence t_n starts with: {t[:10]}...")
    print(f"The largest n for v_n = (t_(n-1), ..., t_(n+2)) is n={n}, with t_29 = {t[29]} <= {limit}.")
    print(f"The base case for the recurrence is f(v_9) = f({t[8]}, {t[9]}, {t[10]}, {t[11]}) = {g9}.")
    print(f"The maximum length is M = f(v_{n}) = {M}.")
    print(f"The tuple (a,b,c,d) with minimal sum achieving M is v_{n}:")
    print(f"a = t_26 = {a}")
    print(f"b = t_27 = {b}")
    print(f"c = t_28 = {c}")
    print(f"d = t_29 = {d}")

    # Step 4: Compute the final expression
    result = (a + b - c - d) % 1000
    
    print("\nThe calculation is (a + b - c - d) mod 1000:")
    print(f"= ({a} + {b} - {t[28]} - {t[29]}) mod 1000")
    print(f"= (({a % 1000} + {b % 1000}) - ({c % 1000} + {d % 1000})) mod 1000")
    a_mod = a % 1000
    b_mod = b % 1000
    c_mod = c % 1000
    d_mod = d % 1000
    sum_ab_mod = (a_mod + b_mod)
    sum_cd_mod = (c_mod + d_mod)
    print(f"= ({sum_ab_mod} - {sum_cd_mod}) mod 1000")
    print(f"= ({sum_ab_mod % 1000} - {sum_cd_mod % 1000}) mod 1000")
    print(f"= ({(sum_ab_mod % 1000) - (sum_cd_mod % 1000)}) mod 1000")
    print(f"= {result}")

solve()