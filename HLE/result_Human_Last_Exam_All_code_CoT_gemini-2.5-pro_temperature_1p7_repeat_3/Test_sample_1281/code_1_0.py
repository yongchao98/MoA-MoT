def calculate_maximum_charge():
    """
    Calculates the maximum amount a consumer can be charged based on an estimate
    under Ontario's Consumer Protection Act (10% rule).
    """
    
    # Original written estimate given by the roofer
    original_estimate = 3000.00
    
    # The maximum percentage increase allowed over the estimate
    allowed_percentage_increase = 0.10
    
    # Calculate the dollar value of the maximum allowed increase
    increase_amount = original_estimate * allowed_percentage_increase
    
    # Calculate the maximum legal invoice amount
    maximum_legal_charge = original_estimate + increase_amount
    
    # Print the equation as requested
    print(f"Under the Ontario Consumer Protection Act, a price cannot exceed an estimate by more than 10%.")
    print(f"Calculation for Marc's maximum payable amount:")
    # The final output needs to be the equation itself, as per the instructions
    print(f"{original_estimate:.2f} + ({original_estimate:.2f} * {allowed_percentage_increase:.2f}) = {maximum_legal_charge:.2f}")

calculate_maximum_charge()