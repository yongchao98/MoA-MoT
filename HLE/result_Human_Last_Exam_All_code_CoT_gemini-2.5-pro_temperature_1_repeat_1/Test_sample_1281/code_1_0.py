def calculate_consumer_obligation():
    """
    Calculates the amount a consumer is required to pay based on an estimate
    under the Ontario Consumer Protection Act's 10% rule.
    """
    estimate_price = 3000.00
    invoiced_price = 3500.00
    allowed_overage_percent = 0.10

    # Calculate the dollar value of the allowed overage
    overage_amount = estimate_price * allowed_overage_percent

    # Calculate the maximum price that can be legally charged
    max_legal_price = estimate_price + overage_amount

    print("Analyzing Marc's payment obligation based on Ontario's Consumer Protection Act.")
    print(f"Original Written Estimate: ${estimate_price:.2f}")
    print(f"Invoiced Amount: ${invoiced_price:.2f}")
    print("\nThe law states a supplier cannot charge more than 10% above their estimate.")
    print("-" * 60)
    print("Step 1: Calculate the maximum allowable overage amount.")
    print(f"Calculation: ${estimate_price:.2f} (estimate) * {allowed_overage_percent:.2f} (10%) = ${overage_amount:.2f}")
    print("\nStep 2: Calculate the total maximum legally chargeable price.")
    print(f"Calculation: ${estimate_price:.2f} (estimate) + ${overage_amount:.2f} (10% overage) = ${max_legal_price:.2f}")
    print("-" * 60)
    print(f"\nConclusion: Since the invoice for ${invoiced_price:.2f} exceeds the maximum legal price of ${max_legal_price:.2f},")
    print(f"Marc is only required to pay HR ${max_legal_price:.2f}.")

calculate_consumer_obligation()