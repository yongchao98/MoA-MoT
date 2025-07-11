def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.
    
    Variables:
    q: initial number of outstanding shares
    E: total market value of equity at time 0
    d: total dividends in year 1 under the original policy
    g: annual dividend growth rate (as a decimal)
    """

    # --- User-defined variables ---
    # You can change these values to match a specific scenario.
    q = 1000000  # Number of shares
    E = 50000000  # Total market value of equity ($)
    d = 2000000   # Original total dividends in year 1 ($)
    g = 0.05      # Annual dividend growth rate (5%)

    # The formula for the per-share ex-dividend price p_1 in year 1 is:
    # p_1 = (E * (1 + g) - d) / q
    
    # Calculate the per-share ex-dividend price
    p_1 = (E * (1 + g) - d) / q

    # Print the explanation and the result
    print("The formula for the per-share ex-dividend price (p_1) in year 1 is:")
    print("p_1 = (E * (1 + g) - d) / q\n")
    
    print("Substituting the given values:")
    print(f"p_1 = ({E:,.0f} * (1 + {g}) - {d:,.0f}) / {q:,.0f}")
    
    intermediate_value = E * (1 + g)
    print(f"p_1 = ({intermediate_value:,.0f} - {d:,.0f}) / {q:,.0f}")
    
    numerator = intermediate_value - d
    print(f"p_1 = {numerator:,.0f} / {q:,.0f}\n")
    
    print(f"The per-share ex-dividend price p_1 is: ${p_1:.2f}")

# Execute the function
calculate_ex_dividend_price()
