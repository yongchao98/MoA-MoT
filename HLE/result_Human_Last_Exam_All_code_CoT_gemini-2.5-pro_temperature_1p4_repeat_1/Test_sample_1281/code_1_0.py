def calculate_roofing_cost():
    """
    Calculates the maximum legal amount Marc has to pay based on the
    Ontario Consumer Protection Act.
    """

    # --- Given Information ---
    estimate_amount = 3000.00
    invoiced_amount = 3500.00
    
    # Under the Ontario Consumer Protection Act, the price cannot exceed the
    # estimate by more than 10%.
    allowable_overage_rate = 0.10

    # --- Plan & Explanation ---
    print("The dispute is governed by Ontario's Consumer Protection Act because the business (HR) is located in Ontario and the contract was made there.")
    print("The Act states that a final invoice cannot exceed a written estimate by more than 10%.")
    print("\nLet's calculate the maximum amount Marc is required to pay.")
    print("-" * 50)
    
    # --- Calculation ---
    allowable_overage_amount = estimate_amount * allowable_overage_rate
    max_legal_price = estimate_amount + allowable_overage_amount

    # --- Output Final Equation and Result ---
    print(f"The original estimate was: ${estimate_amount:.2f}")
    print(f"The allowable overage is 10%, which is: ${allowable_overage_amount:.2f}")
    
    print("\nThe final equation for the maximum legal price is:")
    # Here we output each number in the final equation as requested.
    print(f"Maximum Legal Price = ${estimate_amount:.2f} (Estimate) + ${allowable_overage_amount:.2f} (10% Overage)")
    
    print(f"\nTherefore, the maximum Marc is required to pay is: ${max_legal_price:.2f}")
    print("-" * 50)

calculate_roofing_cost()
<<<B>>>