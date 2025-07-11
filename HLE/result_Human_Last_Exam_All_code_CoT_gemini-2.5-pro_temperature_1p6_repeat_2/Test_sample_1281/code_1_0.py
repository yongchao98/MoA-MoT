def calculate_maximum_charge():
    """
    Calculates the maximum legal charge based on an estimate under
    Ontario's Consumer Protection Act.
    """
    # The written estimate given to Marc.
    estimate_price = 3000.00

    # The final invoiced amount.
    invoiced_price = 3500.00

    # The legally allowed percentage increase over the estimate.
    allowed_overage_percent = 10.0

    # 1. Calculate the dollar value of the allowed overage.
    allowed_overage_amount = estimate_price * (allowed_overage_percent / 100.0)

    # 2. Calculate the maximum legally permissible price.
    maximum_legal_price = estimate_price + allowed_overage_amount

    # 3. Print the breakdown of the calculation as an equation.
    print("Under Ontario's Consumer Protection Act, an invoice cannot exceed a written estimate by more than 10%.")
    print("\nHere is the calculation of the maximum legal price:")
    print("-------------------------------------------------------")
    # This line prints each number in the final equation.
    print(f"Maximum Price = (Original Estimate) + (Original Estimate * 10%)")
    print(f"Maximum Price = ${estimate_price:.2f} + (${estimate_price:.2f} * {allowed_overage_percent/100})")
    print(f"Maximum Price = ${estimate_price:.2f} + ${allowed_overage_amount:.2f}")
    print(f"Maximum Legal Price = ${maximum_legal_price:.2f}")
    print("-------------------------------------------------------")

    print(f"\nThe roofer invoiced Marc for ${invoiced_price:.2f}.")
    print(f"Because this amount is greater than the maximum legal price of ${maximum_legal_price:.2f}, Marc is only required to pay ${maximum_legal_price:.2f}.")

calculate_maximum_charge()
<<<B>>>