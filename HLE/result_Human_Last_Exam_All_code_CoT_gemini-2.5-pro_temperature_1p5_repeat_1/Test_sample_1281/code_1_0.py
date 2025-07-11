def calculate_maximum_charge():
    """
    Calculates the maximum legal price based on an estimate under the Ontario CPA.
    """
    # The written estimate given to Marc
    estimate_amount = 3000.00

    # The percentage by which a price can legally exceed an estimate
    allowed_overage_percent = 0.10

    # Calculate the dollar value of the allowed overage
    overage_amount = estimate_amount * allowed_overage_percent

    # Calculate the total maximum legal price
    maximum_legal_price = estimate_amount + overage_amount

    print("The Ontario Consumer Protection Act limits price increases over an estimate to 10%.")
    print("\nHere is the calculation for Marc's situation:")
    # The final print statement shows each number in the final equation.
    print(f"The maximum legal price is the original estimate of ${estimate_amount:,.2f} plus the 10% allowed overage of ${overage_amount:,.2f}, which equals ${maximum_legal_price:,.2f}.")
    print(f"\nFinal Equation: ${estimate_amount:,.2f} + ${overage_amount:,.2f} = ${maximum_legal_price:,.2f}")


calculate_maximum_charge()
