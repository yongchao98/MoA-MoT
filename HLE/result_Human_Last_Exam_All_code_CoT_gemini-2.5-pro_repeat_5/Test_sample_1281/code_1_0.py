def analyze_roofing_dispute():
    """
    Analyzes Marc's roofing dispute based on the Ontario Consumer Protection Act.
    """
    # --- Case Variables ---
    estimate_price = 3000.00
    invoiced_price = 3500.00
    cpa_limit_percent = 10.0

    # --- Legal Analysis and Calculation ---
    print("Legal Analysis of the Consumer Dispute:")
    print("-" * 40)
    print("The transaction is governed by the Ontario Consumer Protection Act (CPA) because the supplier is an Ontario company.")
    print("Section 10 of the CPA states that a supplier cannot charge more than 10% over a written estimate.")
    print("\nCalculating the maximum legal charge:")
    
    # Calculate the maximum allowable overage amount
    max_overage_amount = estimate_price * (cpa_limit_percent / 100.0)
    
    # Calculate the total maximum legal price
    max_legal_price = estimate_price + max_overage_amount
    
    # Print the equation with all the numbers as requested
    print(f"\nOriginal Estimate + (Estimate * 10%) = Maximum Legal Price")
    print(f"${estimate_price:,.2f} + (${estimate_price:,.2f} * 0.10) = ${max_legal_price:,.2f}")

    # --- Conclusion ---
    print("\nConclusion:")
    print(f"The invoiced price of ${invoiced_price:,.2f} is higher than the maximum legally allowed price of ${max_legal_price:,.2f}.")
    print("Therefore, Marc is only required to pay HR the maximum amount permitted by the CPA.")
    print(f"\nFinal amount Marc must pay: ${max_legal_price:,.2f}")

analyze_roofing_dispute()