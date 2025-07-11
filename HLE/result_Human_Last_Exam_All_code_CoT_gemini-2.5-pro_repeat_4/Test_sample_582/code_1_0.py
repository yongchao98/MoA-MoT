def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the old policy.
        g (float): Dividend growth rate.
    """
    if q <= 0:
        print("Number of shares (q) must be positive.")
        return

    # The derived formula for the per-share ex-dividend price in year 1 is:
    # p1 = (E * (1 + g) - d) / q

    p1 = (E * (1 + g) - d) / q

    # Output the explanation and the calculation
    print("The formula for the per-share ex-dividend price (p1) in year 1 is:")
    print("p1 = (E * (1 + g) - d) / q\n")
    print("Substituting the given values:")
    print(f"p1 = ({E:,.0f} * (1 + {g}) - {d:,.0f}) / {q:,.0f}")
    print(f"p1 = ({E * (1 + g):,.0f} - {d:,.0f}) / {q:,.0f}")
    print(f"p1 = {(E * (1 + g) - d):,.0f} / {q:,.0f}")
    print(f"The per-share ex-dividend price in year 1 is: ${p1:.2f}")


# --- User-defined variables ---
# You can change these values to test with different scenarios.
# q: outstanding shares
q = 1000000
# E: total market value of equity
E = 50000000
# d: total dividends in year 1 (original policy)
d = 2500000
# g: annual dividend growth rate
g = 0.05
# --- End of user-defined variables ---

calculate_ex_dividend_price(q, E, d, g)