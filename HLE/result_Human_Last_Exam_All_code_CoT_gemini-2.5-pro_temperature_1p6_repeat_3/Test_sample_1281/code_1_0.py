def calculate_max_payment():
    """
    Calculates the maximum amount payable based on an estimate
    according to Ontario's Consumer Protection Act.
    """
    # The initial written estimate for the job.
    estimate = 3000.00

    # The final invoiced amount.
    invoice_amount = 3500.00

    # Under the Ontario Consumer Protection Act, the final price cannot
    # exceed the estimate by more than 10%.
    allowed_overage_rate = 0.10

    # Calculate the dollar value of the maximum allowed overage.
    max_overage_in_dollars = estimate * allowed_overage_rate

    # Calculate the maximum legal price that can be charged.
    max_legal_price = estimate + max_overage_in_dollars

    print("The original estimate provided was: $3000")
    print("The final invoice was for: $3500")
    print("\nAccording to the Ontario Consumer Protection Act, the final price cannot exceed the estimate by more than 10%.")
    print(f"The calculation for the maximum legal price is as follows:")
    
    # Printing the final equation with each number.
    print(f"${int(estimate)} (Estimate) + ${int(max_overage_in_dollars)} (10% Overage) = ${int(max_legal_price)} (Max Legal Price)")
    
    print(f"\nSince the invoice for ${int(invoice_amount)} is greater than the maximum legal price of ${int(max_legal_price)}, Marc is only required to pay ${int(max_legal_price)}.")

calculate_max_payment()