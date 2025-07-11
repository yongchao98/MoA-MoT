def calculate_roofing_cost():
    """
    Calculates the maximum legal amount Marc has to pay based on the
    Ontario Consumer Protection Act.
    """
    # Given values from the scenario
    estimate = 3000.00
    invoiced_amount = 3500.00
    allowable_overage_percentage = 0.10  # 10%

    # Explanation of the relevant law
    print("According to Ontario's Consumer Protection Act, a business cannot charge more than 10% above a written estimate.")
    print("-" * 40)
    
    # Step 1: Calculate the maximum allowable overage in dollars
    max_overage_amount = estimate * allowable_overage_percentage
    
    # Step 2: Calculate the maximum legal price
    max_legal_price = estimate + max_overage_amount

    # Step 3: Output the calculation step-by-step, showing all numbers
    print(f"Original Estimate: ${estimate:,.2f}")
    print(f"Allowable Overage: {allowable_overage_percentage * 100}%")
    print("\nFinal Equation:")
    print(f"Maximum Legal Price = Original Estimate + (Original Estimate * Allowable Overage)")
    print(f"Maximum Legal Price = ${estimate:,.2f} + (${estimate:,.2f} * {allowable_overage_percentage})")
    print(f"Maximum Legal Price = ${estimate:,.2f} + ${max_overage_amount:,.2f} = ${max_legal_price:,.2f}")
    print("-" * 40)
    
    # Final conclusion
    print(f"The invoiced amount was ${invoiced_amount:,.2f}, which is more than the maximum legal price.")
    print(f"Therefore, Marc is only legally required to pay HR ${max_legal_price:,.2f}.")

# Run the calculation
calculate_roofing_cost()
<<<B>>>