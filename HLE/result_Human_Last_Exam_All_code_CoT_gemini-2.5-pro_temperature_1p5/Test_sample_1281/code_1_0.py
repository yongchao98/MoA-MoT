def calculate_final_price():
    """
    Calculates the maximum legal price based on Ontario's Consumer Protection Act
    and explains the reasoning.
    """
    estimate = 3000.00
    invoice = 3500.00
    overage_limit_percentage = 0.10  # 10%

    # Calculate the maximum legally allowed overage
    max_overage_amount = estimate * overage_limit_percentage

    # Calculate the total maximum legal price
    max_legal_price = estimate + max_overage_amount

    print("The legal issue concerns Ontario's Consumer Protection Act (CPA) and its rules on estimates.")
    print(f"The original estimate was ${int(estimate)}.")
    print("The CPA prohibits a supplier from charging more than 10% above their estimate.")
    print("\nHere is the calculation for the maximum legally permissible price:")
    
    # As requested, printing each number in the final equation
    print(f"{int(estimate)} + ({int(estimate)} * {overage_limit_percentage}) = {int(max_legal_price)}")
    
    print(f"\nAlthough the final invoice was for ${int(invoice)}, Marc is only required to pay ${int(max_legal_price)}.")
    print("This outcome aligns with the provisions of the Consumer Protection Act.")

calculate_final_price()