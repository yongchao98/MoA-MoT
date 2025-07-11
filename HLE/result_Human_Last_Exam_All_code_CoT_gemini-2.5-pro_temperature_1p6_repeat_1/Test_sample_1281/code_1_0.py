def calculate_maximum_chargeable_amount():
    """
    Calculates the maximum legal price based on an estimate under
    Ontario's Consumer Protection Act.
    """
    # The original written estimate provided by Honest Roofers Inc.
    estimate_price = 3000.00

    # The final amount invoiced to Marc.
    invoiced_price = 3500.00

    # Under the CPA, the final price cannot exceed the estimate by more than 10%.
    allowable_overage_percentage = 0.10

    # Step 1: Calculate the dollar value of the 10% allowable overage.
    max_overage_amount = estimate_price * allowable_overage_percentage

    # Step 2: Calculate the maximum legal price by adding the overage to the estimate.
    max_legal_price = estimate_price + max_overage_amount

    print("Ontario's Consumer Protection Act limits how much a final price can exceed a written estimate.")
    print("The rule is that the final charge cannot be more than 10% above the estimate.")
    print("-" * 50)
    print("Here is the calculation based on Marc's situation:")
    print(f"Original Estimate: ${estimate_price:.2f}")
    print(f"Invoiced Amount: ${invoiced_price:.2f}")

    print("\nCalculating the maximum legal price:")
    # This line prints the final equation with each number, as requested.
    print(f"Equation: ${estimate_price:.2f} (Estimate) + ${max_overage_amount:.2f} (10% Overage) = ${max_legal_price:.2f} (Maximum Legal Price)")

    print("-" * 50)
    print(f"The final invoice of ${invoiced_price:.2f} exceeds the maximum legal price of ${max_legal_price:.2f}.")
    print(f"Therefore, the lawyer would advise Marc that he is only required to pay HR ${max_legal_price:.2f}.")
    print("This corresponds to answer choice B.")

# Execute the function to get the advice for Marc.
calculate_maximum_chargeable_amount()