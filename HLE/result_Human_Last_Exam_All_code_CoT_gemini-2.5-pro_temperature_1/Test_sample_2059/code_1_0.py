def solve():
    """
    Calculates the sum of the squares of the coefficients of the given polynomial product.
    Let S_n be the sum of squares of coefficients for the product up to i=n.
    Let C_n be the sum of products of adjacent coefficients.
    We derived the recurrence relations:
    S_n = 4*S_{n-1} + 2*C_{n-1}
    C_n = 3*S_{n-1} + 2*C_{n-1}
    This leads to a single recurrence for S_n: S_n = 6*S_{n-1} - 2*S_{n-2}.
    The solution can be expressed as S_n = 2*u_n + 5*v_n, where u_n and v_n
    are integer sequences also following x_n = 6*x_{n-1} - 2*x_{n-2}.
    """
    n = 19

    # Initial values for the u_n sequence
    # u_n = (3+sqrt(7))^n + (3-sqrt(7))^n
    u_prev = 2  # u_0
    u_curr = 6  # u_1

    # Initial values for the v_n sequence
    # v_n = ((3+sqrt(7))^n - (3-sqrt(7))^n) / sqrt(7)
    v_prev = 0  # v_0
    v_curr = 2  # v_1

    if n == 0:
        u_n = u_prev
        v_n = v_prev
    elif n == 1:
        u_n = u_curr
        v_n = v_curr
    else:
        for _ in range(2, n + 1):
            u_new = 6 * u_curr - 2 * u_prev
            u_prev = u_curr
            u_curr = u_new

            v_new = 6 * v_curr - 2 * v_prev
            v_prev = v_curr
            v_curr = v_new
        u_n = u_curr
        v_n = v_curr

    # The final sum S_19 is given by 2*u_19 + 5*v_19
    result = 2 * u_n + 5 * v_n

    print(f"The calculation is based on the formula: 2 * u_{n} + 5 * v_{n}")
    print(f"u_{n} = {u_n}")
    print(f"v_{n} = {v_n}")
    print(f"Result = 2 * {u_n} + 5 * {v_n} = {result}")

solve()