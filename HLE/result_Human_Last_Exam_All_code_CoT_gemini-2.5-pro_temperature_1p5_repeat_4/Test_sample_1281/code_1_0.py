def calculate_maximum_charge():
    """
    Calculates the maximum legal charge based on an estimate
    under Ontario's Consumer Protection Act.
    """
    # Original written estimate given to Marc
    estimated_price = 3000.00

    # Invoiced price after the job was completed
    invoiced_price = 3500.00

    # The maximum allowable increase over the estimate is 10%
    allowable_increase_factor = 1.10

    # Calculate the maximum price HR can legally charge
    maximum_allowable_price = estimated_price * allowable_increase_factor

    # Explain the legal rule and show the calculation
    print("Under Ontario's Consumer Protection Act, the final price cannot exceed the estimate by more than 10%.")
    print("\nHere is the calculation based on the provided numbers:")
    print(f"The original estimate was ${int(estimated_price):,}.00.")
    print(f"The maximum allowable increase is 10%.")
    print(f"\nThe equation for the maximum legal price is:")
    print(f"${int(estimated_price):,}.00 * 1.10 = ${int(maximum_allowable_price):,}.00")
    print(f"\nSince the invoiced price of ${int(invoiced_price):,}.00 exceeds this amount, the lawyer would advise Marc that he is only required to pay ${int(maximum_allowable_price):,}.00.")

calculate_maximum_charge()