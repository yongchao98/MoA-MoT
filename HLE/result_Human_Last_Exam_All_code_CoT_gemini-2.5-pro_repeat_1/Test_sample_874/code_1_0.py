def solve():
    """
    Finds the specific tuple (a, b, c, d) and computes the required expression.
    """
    
    # Step 1: Define the Tribonacci sequence t_n starting t_0=0, t_1=0, t_2=1.
    # We need to find the largest n such that t_{n+2} <= 10,000,000.
    
    limit = 10_000_000
    
    t = [0, 0, 1]
    while t[-1] <= limit:
        t.append(t[-1] + t[-2] + t[-3])
    
    # The last element of t is the first one to exceed the limit.
    # So, the largest Tribonacci number we can use is t[-2].
    # This corresponds to t_{n+2}. The index of t[-2] is len(t)-2.
    # So, n+2 = len(t)-2, which means n = len(t)-4.
    
    n = len(t) - 4
    
    # The tuple is (a,b,c,d) = (t_{n-1}-1, t_n, t_{n+1}, t_{n+2})
    # which corresponds to indices n-1, n, n+1, n+2 in the generated list `t`.
    # Let's verify the indices
    # k = n + 2
    # t_k = t[k]
    # a = t[k-3]-1, b=t[k-2], c=t[k-1], d=t[k]
    
    k = n + 2 # index for d

    # Step 2: Compute a, b, c, d modulo 1000.
    # We can do this by generating the Tribonacci sequence modulo 1000.
    m = 1000
    t_mod = [0, 0, 1]
    # We need up to index k.
    for i in range(3, k + 1):
        next_t = (t_mod[i-1] + t_mod[i-2] + t_mod[i-3]) % m
        t_mod.append(next_t)
        
    a_mod = (t_mod[k-3] - 1 + m) % m
    b_mod = t_mod[k-2]
    c_mod = t_mod[k-1]
    d_mod = t_mod[k]

    # For verification, let's print the actual values.
    # print(f"Found n = {n}")
    # print(f"Using Tribonacci numbers with indices {k-3, k-2, k-1, k}")
    a_val = t[k-3] - 1
    b_val = t[k-2]
    c_val = t[k-1]
    d_val = t[k]
    # print(f"Tuple (a,b,c,d) = ({a_val}, {b_val}, {c_val}, {d_val})")

    # Step 3: Compute the expression (a + b - c - d) mod 1000
    result = (a_mod + b_mod - c_mod - d_mod + m) % m
    
    # Output the final equation with each number
    print(f"The determined tuple (a, b, c, d) modulo 1000 is ({a_mod}, {b_mod}, {c_mod}, {d_mod}).")
    print(f"The expression to compute is (a + b - c - d) mod 1000.")
    print(f"Calculation: ({a_mod} + {b_mod} - {c_mod} - {d_mod}) mod 1000 = {result}")
    print(f"The final answer is: {result}")

solve()