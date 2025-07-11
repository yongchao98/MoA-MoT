def calculate_maximum_charge():
    """
    Calculates the maximum legal charge based on a written estimate
    under Ontario's Consumer Protection Act.
    """
    # The written estimate provided to Marc.
    original_estimate = 3000.00

    # The percentage by which a final price can legally exceed an estimate.
    allowable_increase_percentage = 0.10  # 10%

    # Calculate the dollar value of the allowable increase.
    allowable_increase_amount = original_estimate * allowable_increase_percentage

    # Calculate the maximum price that can be legally charged.
    max_legal_price = original_estimate + allowable_increase_amount

    print("According to Ontario's Consumer Protection Act, the final price cannot exceed the written estimate by more than 10%.")
    print("Here is the calculation for Marc's situation:")
    print("-" * 40)
    print(f"Original Written Estimate: ${original_estimate:.2f}")
    print(f"Allowable 10% Increase Amount: ${allowable_increase_amount:.2f}")
    print("-" * 40)
    print("The maximum legal price is calculated as follows:")
    # Final equation with all numbers printed out
    print(f"${original_estimate:.2f} (Estimate) + ${allowable_increase_amount:.2f} (10% Buffer) = ${max_legal_price:.2f} (Maximum Legal Price)")
    print("\nTherefore, the lawyer would advise Marc that he is only required to pay $3,300.00.")

calculate_maximum_charge()
<<<B>>>