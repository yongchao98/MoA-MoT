def calculate_roofing_cost():
    """
    Calculates the maximum legal price Marc must pay based on the Ontario Consumer Protection Act.
    """
    
    # Initial values from the scenario
    original_estimate = 3000.00
    invoiced_price = 3500.00
    
    # The Ontario Consumer Protection Act allows for a maximum 10% overage on estimates
    cpa_allowed_overage_percent = 0.10
    
    # Calculate the maximum dollar amount of the overage
    max_overage_amount = original_estimate * cpa_allowed_overage_percent
    
    # Calculate the total maximum legal price
    max_legal_price = original_estimate + max_overage_amount
    
    print("The lawyer would advise that Ontario's Consumer Protection Act applies.")
    print("This Act states that a final price cannot exceed a written estimate by more than 10%.")
    print("\nHere is the calculation for the maximum amount Marc is required to pay:")
    
    # The final print statement shows each number in the final equation.
    print(f"Original Estimate + (10% of Estimate) = Maximum Legal Price")
    print(f"${original_estimate:.2f} + ${max_overage_amount:.2f} = ${max_legal_price:.2f}")

    print(f"\nTherefore, Marc is only required to pay HR $3,300.00, as stated in answer choice B.")

calculate_roofing_cost()
<<<B>>>