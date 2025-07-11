import math

def get_odd_part(m):
    """
    Calculates the odd part of a positive integer m.
    This is done by repeatedly dividing m by 2 until it is odd.
    """
    if m == 0:
        return 0
    while m % 2 == 0:
        m //= 2
    return m

def S(m):
    """
    Calculates the number of weighings for the subproblem with m items.
    S(m) = 0 if m is odd.
    S(m) = m - odd_part(m) if m is even.
    """
    if m % 2 != 0:
        return 0
    return m - get_odd_part(m)

def T(n):
    """
    Calculates the minimum number of trials T(n).
    T(n) = n + max(S(m)) for 0 <= m <= n.
    """
    max_s = 0
    # We only need to check even m, as S(m) is 0 for odd m.
    for m in range(0, n + 1, 2):
        s_m = S(m)
        if s_m > max_s:
            max_s = s_m
    return n + max_s

def solve():
    """
    Calculates and prints the values for T(2), T(3), T(1234), and T(6712).
    """
    n_values = [2, 3, 1234, 6712]
    results = [T(n) for n in n_values]
    print(','.join(map(str, results)))

solve()