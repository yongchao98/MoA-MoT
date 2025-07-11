def calculate_roofing_cost():
    """
    Calculates the maximum legal price for a service based on an estimate,
    according to Ontario's Consumer Protection Act.
    """
    # Original details of the transaction
    original_estimate = 3000.00
    invoiced_price = 3500.00
    
    # Under the Ontario Consumer Protection Act, the final price cannot exceed
    # the estimate by more than 10%.
    allowable_increase_percentage = 10.0
    
    # Calculate the maximum allowable price
    allowable_multiplier = 1 + (allowable_increase_percentage / 100.0)
    max_legal_price = original_estimate * allowable_multiplier
    
    print("This scenario is governed by Ontario's Consumer Protection Act.")
    print("The Act states that a final invoice cannot exceed the written estimate by more than 10%.")
    print("-" * 60)
    
    # Show the final equation with all the numbers
    print("Calculation of the maximum allowable charge:")
    print(f"Original Estimate: ${original_estimate:,.2f}")
    print(f"Maximum Increase Percentage: {allowable_increase_percentage}%")
    print(f"Final Equation: ${original_estimate:,.2f} * (1 + {allowable_increase_percentage / 100}) = ${max_legal_price:,.2f}")
    print("-" * 60)
    
    # State the conclusion
    print(f"Therefore, although the invoice was for ${invoiced_price:,.2f}, Marc is only legally required to pay ${max_legal_price:,.2f}.")

calculate_roofing_cost()
<<<B>>>