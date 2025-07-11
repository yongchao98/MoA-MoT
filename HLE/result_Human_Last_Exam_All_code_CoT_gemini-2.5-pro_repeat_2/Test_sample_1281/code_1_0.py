def calculate_maximum_charge():
    """
    Calculates the maximum legal price based on a written estimate and the 
    Ontario Consumer Protection Act's 10% rule.
    """
    # Original written estimate given by Honest Roofers Inc.
    estimated_price = 3000.00

    # The maximum percentage increase allowed over an estimate by the Consumer Protection Act.
    allowable_increase_percentage = 0.10

    # Calculate the dollar value of the maximum allowable increase.
    max_increase_amount = estimated_price * allowable_increase_percentage

    # Calculate the total maximum legal price.
    max_legal_price = estimated_price + max_increase_amount

    print("Under Ontario's Consumer Protection Act, the final price cannot exceed a written estimate by more than 10%.")
    print("\nCalculation of the maximum legal price:")
    print(f"Original Estimate: ${estimated_price:.2f}")
    print(f"Maximum Allowable Increase: ${estimated_price:.2f} * {allowable_increase_percentage:.2f} = ${max_increase_amount:.2f}")
    print(f"Maximum Legal Price Marc is required to pay: ${estimated_price:.2f} + ${max_increase_amount:.2f} = ${max_legal_price:.2f}")
    print("\nThis means the correct course of action aligns with option B.")

calculate_maximum_charge()