def calculate_final_bill():
    """
    Calculates the maximum payable amount based on an estimate and a percentage limit.
    """
    original_estimate = 3000.00
    invoiced_amount = 3500.00
    allowed_increase_percentage = 0.10  # 10%

    # Calculate the maximum allowed increase in dollars
    max_increase_amount = original_estimate * allowed_increase_percentage

    # Calculate the total legally payable amount
    total_payable = original_estimate + max_increase_amount

    # Print the explanation and the final equation
    print("According to Ontario's Consumer Protection Act, a final invoice cannot exceed the written estimate by more than 10%.")
    print("\nCalculation:")
    print(f"Original Estimate: ${original_estimate:.2f}")
    print(f"Maximum Increase (10%): ${max_increase_amount:.2f}")
    print("The final amount Marc is required to pay is calculated as follows:")
    print(f"${original_estimate:.2f} (Estimate) + ${max_increase_amount:.2f} (10% Increase) = ${total_payable:.2f}")

calculate_final_bill()