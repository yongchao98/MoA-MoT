def calculate_maximum_charge():
    """
    Calculates the maximum legal price based on an estimate and the
    Ontario Consumer Protection Act's 10% rule.
    """
    # Original written estimate given to Marc
    original_estimate = 3000.00

    # The maximum allowable percentage increase over an estimate
    allowable_percentage_increase = 0.10

    # Step 1: Calculate the dollar value of the allowable 10% increase
    allowable_overage_amount = original_estimate * allowable_percentage_increase

    # Step 2: Calculate the maximum legal price
    maximum_legal_price = original_estimate + allowable_overage_amount

    print("Analyzing the charge based on the Ontario Consumer Protection Act:")
    print("-" * 60)
    print(f"Original Written Estimate: ${original_estimate:.2f}")
    print(f"Invoiced Amount: $3500.00")
    print(f"Allowable Overage Percentage: {allowable_percentage_increase * 100:.0f}%")
    print("-" * 60)
    print("The maximum legal price is calculated by adding the 10% allowable overage to the original estimate.")
    print("\nHere is the final equation showing each number:")
    print(f"${original_estimate:.2f} (Original Estimate) + ${allowable_overage_amount:.2f} (10% Overage) = ${maximum_legal_price:.2f} (Maximum Legal Price)")
    print(f"\nConclusion: Marc is required to pay HR ${maximum_legal_price:.2f}.")

# Run the calculation and print the result
calculate_maximum_charge()