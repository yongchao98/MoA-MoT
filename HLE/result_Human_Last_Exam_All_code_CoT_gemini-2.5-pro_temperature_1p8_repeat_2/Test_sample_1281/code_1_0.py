def calculate_maximum_charge():
    """
    Calculates the maximum legal charge based on an estimate under Ontario's
    Consumer Protection Act, which limits price increases to 10% over the estimate.
    """
    original_estimate = 3000.00
    percent_limit = 0.10
    
    # Calculate the maximum allowable increase
    allowed_increase_amount = original_estimate * percent_limit
    
    # Calculate the maximum total legal charge
    maximum_legal_charge = original_estimate + allowed_increase_amount
    
    # The actual invoice presented to Marc
    invoiced_amount = 3500.00

    print(f"The original written estimate was: ${original_estimate:,.2f}")
    print(f"Under the Consumer Protection Act, the price cannot exceed the estimate by more than 10%.")
    print(f"The maximum allowed increase is 10% of ${original_estimate:,.2f}, which is ${allowed_increase_amount:,.2f}.")
    print(f"Therefore, the maximum legal amount Marc is required to pay is: ${original_estimate:,.2f} + ${allowed_increase_amount:,.2f} = ${maximum_legal_charge:,.2f}.")
    print(f"The invoice for ${invoiced_amount:,.2f} is higher than the legally permissible amount.")

calculate_maximum_charge()