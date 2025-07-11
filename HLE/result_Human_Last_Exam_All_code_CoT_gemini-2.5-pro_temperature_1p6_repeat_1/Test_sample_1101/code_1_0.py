import math

def solve():
    """
    Calculates the sum S_n and demonstrates its bound f(n) = n! * 2^n.
    S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(3/2+n)/Gamma(3/2+n-m)
    """

    max_n = 20
    s = {0: 1.0, 1: -0.5}

    print("n, S_n, f(n), |S_n/f(n)|")
    print("-" * 35)

    # Calculate for n=0
    f_0 = math.factorial(0) * (2**0)
    ratio_0 = abs(s[0] / f_0)
    print(f"{0:2d}, {s[0]:10.4f}, {f_0:10.1f}, {ratio_0:10.6f}")

    # Calculate for n=1
    f_1 = math.factorial(1) * (2**1)
    ratio_1 = abs(s[1] / f_1)
    print(f"{1:2d}, {s[1]:10.4f}, {f_1:10.1f}, {ratio_1:10.6f}")

    # Calculate for n > 1 using the recurrence relation
    # 2*S_{n+1} + (4n+1)*S_n + n*(2n+1)*S_{n-1} = 0
    # S_{n+1} = -((4n+1)*S_n + n*(2n+1)*S_{n-1}) / 2
    for n in range(1, max_n):
        s_n_plus_1 = -((4 * n + 1) * s[n] + n * (2 * n + 1) * s[n - 1]) / 2.0
        s[n + 1] = s_n_plus_1

        f_n_plus_1 = math.factorial(n + 1) * (2**(n + 1))
        ratio = abs(s[n + 1] / f_n_plus_1)
        
        # Output each number in the final inequality relation demonstration
        # |S_n| <= C * f(n)  or  |S_n / f(n)| <= C
        # Here we show the value of |S_n / f(n)| which is bounded
        print(f"{n+1:2d}, {s[n+1]:10.4f}, {f_n_plus_1:10.1e}, {ratio:10.6f}")

solve()