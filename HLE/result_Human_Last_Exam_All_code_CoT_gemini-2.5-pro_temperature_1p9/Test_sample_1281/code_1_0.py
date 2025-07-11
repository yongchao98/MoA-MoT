def calculate_max_payment():
    """
    Calculates the maximum legal amount payable based on the Ontario Consumer Protection Act's
    rule on estimates.
    """
    initial_estimate = 3000.00
    allowed_percentage_increase = 0.10

    # Calculate the maximum allowed increase in dollar value
    max_increase_amount = initial_estimate * allowed_percentage_increase

    # Calculate the total maximum legal price
    max_legal_price = initial_estimate + max_increase_amount

    print("Analyzing the payment dispute based on the Ontario Consumer Protection Act:")
    print("-" * 70)
    print(f"Original written estimate: ${initial_estimate:,.2f}")
    print(f"The Act states that the final price cannot exceed the estimate by more than {allowed_percentage_increase:.0%}.")
    print("\nCalculation of the maximum legal charge:")
    # The final equation as requested, showing each number
    print(f"Maximum Legal Price = Original Estimate + (Original Estimate * {allowed_percentage_increase:.0%})")
    print(f"Maximum Legal Price = ${initial_estimate:,.2f} + ${max_increase_amount:,.2f}")
    print(f"Total maximum Marc is required to pay: ${max_legal_price:,.2f}")

calculate_max_payment()