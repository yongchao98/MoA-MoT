def calculate_max_charge():
    """
    Calculates the maximum legal charge based on an estimate and a percentage limit,
    as per Ontario's Consumer Protection Act.
    """
    estimate_price = 3000.00
    max_increase_percentage = 10.0
    
    # Calculate the dollar value of the maximum allowed increase
    max_increase_amount = estimate_price * (max_increase_percentage / 100.0)
    
    # Calculate the total maximum legal price
    max_legal_price = estimate_price + max_increase_amount
    
    print("Based on Ontario's Consumer Protection Act, the final price cannot exceed the estimate by more than 10%.")
    print("\n--- Calculation ---")
    print(f"Original Estimate: ${estimate_price:.2f}")
    print(f"Maximum Allowed Increase: {max_increase_percentage}%")
    print(f"Maximum Increase Amount: ${max_increase_amount:.2f}")
    print("\n--- Final Price Equation ---")
    print(f"${estimate_price:.2f} (Estimate) + ${max_increase_amount:.2f} (10% Increase) = ${max_legal_price:.2f}")
    print("\nTherefore, the maximum amount Marc is legally required to pay is $3,300.00.")

calculate_max_charge()
<<<B>>>