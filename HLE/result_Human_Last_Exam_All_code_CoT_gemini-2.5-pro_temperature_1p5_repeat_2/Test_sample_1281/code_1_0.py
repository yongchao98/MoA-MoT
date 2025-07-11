def calculate_amount_owed():
    """
    Calculates the amount Marc owes based on Ontario's Consumer Protection Act.
    The Act states a final price cannot exceed a written estimate by more than 10%.
    """
    
    # Initial figures from the problem
    estimate_price = 3000.00
    invoiced_price = 3500.00
    allowable_increase_percent = 0.10

    # Step 1: Explain the legal principle
    print("According to Ontario's Consumer Protection Act, the final price cannot exceed a written estimate by more than 10%.")
    
    # Step 2: Calculate the maximum allowed increase in dollars
    max_increase_amount = estimate_price * allowable_increase_percent
    
    # Step 3: Calculate the maximum price that can be legally charged
    max_legal_price = estimate_price + max_increase_amount
    
    # Step 4: Determine the final amount owed. Marc must pay the maximum legal price
    # because the invoiced amount ($3,500) is higher.
    amount_owed = max_legal_price

    print("\nHere is the calculation based on the law:")
    print("------------------------------------------")
    # Step 5: Output the final equation with all its parts, as requested.
    print(f"Original Estimate: ${estimate_price:,.2f}")
    print(f"Maximum 10% Increase: ${max_increase_amount:,.2f}")
    print(f"Invoice from Honest Roofers Inc.: ${invoiced_price:,.2f}")
    print("\nThe final amount Marc is required to pay is calculated as:")
    print(f"${estimate_price:,.2f} (Estimate) + ${max_increase_amount:,.2f} (10% Increase) = ${amount_owed:,.2f} (Amount Owed)")
    print("------------------------------------------")

# Run the calculation
calculate_amount_owed()
<<<B>>>