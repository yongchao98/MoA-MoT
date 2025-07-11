def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the ex-dividend price per share in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the old policy.
        g (float): Annual dividend growth rate.
    """
    # The formula derived is p1 = (E * (1 + g) - d) / q
    p1 = (E * (1 + g) - d) / q

    # Print the equation with the given values
    print("The formula for the ex-dividend price per share (p1) is:")
    print("p1 = (E * (1 + g) - d) / q")
    print("\nSubstituting the given values:")
    # Use format specifiers to avoid potential floating point inaccuracies in the print string
    print(f"p1 = ({E:g} * (1 + {g:g}) - {d:g}) / {q:g}")

    # Calculate and print the final result
    print("\nResult:")
    print(f"p1 = {p1}")

# --- User-provided inputs ---
# Number of outstanding shares
q = 10000000
# Total market value of equity
E = 500000000
# Total dividends in year 1 (original policy)
d = 25000000
# Annual dividend growth rate
g = 0.06
# --- End of inputs ---

calculate_ex_dividend_price(q, E, d, g)