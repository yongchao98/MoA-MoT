def calculate_maximum_chargeable_amount():
    """
    Calculates the maximum amount a consumer can be charged based on an estimate
    under Ontario's Consumer Protection Act and prints the result.
    """
    initial_estimate = 3000.00
    invoiced_amount = 3500.00
    overage_limit_percent = 10

    # Calculate the maximum allowed overage in dollars
    max_overage_decimal = overage_limit_percent / 100.0
    max_overage_amount = initial_estimate * max_overage_decimal
    
    # Calculate the maximum legal price
    max_legal_price = initial_estimate + max_overage_amount

    print("Analyzing the price dispute based on Ontario's Consumer Protection Act:")
    print(f"1. The original written estimate was: ${initial_estimate:,.2f}")
    print(f"2. The Act allows for the final price to exceed the estimate by no more than {overage_limit_percent}%.")
    print("\nCalculating the maximum legally-allowed price:")
    print(f"The maximum allowed overage is {overage_limit_percent}% of ${initial_estimate:,.2f}, which is ${max_overage_amount:,.2f}.")
    print("\nThe final equation for the maximum price is:")
    print(f"Original Estimate + Maximum Overage = Maximum Legal Price")
    print(f"${initial_estimate:,.2f} + ${max_overage_amount:,.2f} = ${max_legal_price:,.2f}")

    print(f"\nConclusion:")
    print(f"The invoiced amount was ${invoiced_amount:,.2f}, which is higher than the allowed maximum.")
    print(f"Therefore, the maximum amount Marc is required to pay is ${max_legal_price:,.2f}.")

calculate_maximum_chargeable_amount()