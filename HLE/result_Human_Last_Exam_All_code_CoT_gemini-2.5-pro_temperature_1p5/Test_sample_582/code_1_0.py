import sys

def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): The number of outstanding shares.
        E (float): The total market value of equity.
        d (float): The total dividends to be distributed in year 1 under the old policy.
        g (float): The annual growth rate of dividends.
    """
    # The derived formula for the per-share ex-dividend price in year 1 is:
    # p_1 = (E * (1 + g) - d) / q
    
    # Check for division by zero
    if q == 0:
        print("Error: The number of shares (q) cannot be zero.", file=sys.stderr)
        return

    # Calculate the price
    p1 = (E * (1 + g) - d) / q
    
    # Print the formula and the calculation steps
    print("The formula for the per-share ex-dividend price (p_1) is:")
    print("p_1 = (E * (1 + g) - d) / q")
    print("\nPlugging in the given values:")
    # Using f-string formatting to display the equation with numbers
    print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p_1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p_1 = {E * (1 + g) - d} / {q}")
    print(f"\nThe per-share ex-dividend price p_1 is: {p1:.2f}")


# Example values from the problem description (using variables)
# You can change these values to see how the price changes.
q_shares = 10000000  # 10 million shares
E_value = 500000000  # $500 million
d_dividends = 20000000  # $20 million
g_growth = 0.05        # 5% growth rate

# Execute the calculation
calculate_ex_dividend_price(q_shares, E_value, d_dividends, g_growth)