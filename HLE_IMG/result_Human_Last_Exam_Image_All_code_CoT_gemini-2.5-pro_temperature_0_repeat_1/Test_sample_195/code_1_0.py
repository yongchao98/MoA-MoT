def solve():
    """
    This function generates and prints the equation for the given graph f(x).
    """
    # Define the factors of the numerator based on the x-intercepts at -b, b, and d.
    # (x - (-b)) -> (x + b)
    # (x - b)
    # (x - d)
    # We can group (x+b) and (x-b) into (x^2 - b^2)
    numerator_str = "(x - d)(x + b)(x - b)"

    # Define the factors of the denominator based on the vertical asymptotes at a and c.
    # (x - a)
    # (x - c)
    denominator_str = "(x - a)(x - c)"

    # Print the final equation for f(x).
    # The equation is constructed to have the lowest possible polynomial order
    # that satisfies the given roots, asymptotes, and end behavior.
    print(f"f(x) = ({numerator_str}) / ({denominator_str})")

solve()
<<<f(x) = ((x - d)(x + b)(x - b)) / ((x - a)(x - c))>>>