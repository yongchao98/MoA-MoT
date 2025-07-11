def solve():
    """
    This function generates and prints the two additional inequalities.
    """
    # The two inequalities are derived using the big-M method to model the piecewise nature
    # of the function f(x) for the case x < 1.
    # The variable 'a' is used to distinguish between x >= 1 and x < 1.
    # The variable 'b' is used to distinguish between x >= 0 and x < 0 when x < 1.
    # A large constant M is required. We choose M=1000 as a typical large value.

    M = 1000

    # First inequality: y >= M*b - M - M*a
    # This corresponds to y >= 0 when a=0, b=1.
    # y >= -M*a + M*b - M
    print(f"y >= {-M}a + {M}b + {-M}")

    # Second inequality: y >= x - M*b - M*a
    # This corresponds to y >= x when a=0, b=0.
    # y >= x - M*a - M*b
    print(f"y >= x + {-M}a + {-M}b")

solve()