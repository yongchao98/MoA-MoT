def calculate_maximum_charge():
    """
    Calculates the maximum legal amount to be charged based on an estimate
    under Ontario's Consumer Protection Act.
    """
    original_estimate = 3000.00
    allowable_overage_percentage = 0.10  # This represents 10%

    # Calculate the overage amount
    overage_amount = original_estimate * allowable_overage_percentage

    # Calculate the maximum legal price
    maximum_legal_price = original_estimate + overage_amount

    # Print the explanation and the final equation
    print("Analyzing the situation based on Ontario's Consumer Protection Act:")
    print("The Act states that the final price cannot exceed the estimate by more than 10%.")
    print("\nHere is the calculation:")
    print("Original Estimate: ${:,.2f}".format(original_estimate))
    print("Allowable Overage Percentage: {}%".format(int(allowable_overage_percentage * 100)))
    print("Final Equation: ${:,.2f} + (${:,.2f} * {}) = ${:,.2f}".format(original_estimate, original_estimate, allowable_overage_percentage, maximum_legal_price))
    print("\nTherefore, the maximum amount Marc is required to pay is ${:,.2f}.".format(maximum_legal_price))

calculate_maximum_charge()
<<<B>>>