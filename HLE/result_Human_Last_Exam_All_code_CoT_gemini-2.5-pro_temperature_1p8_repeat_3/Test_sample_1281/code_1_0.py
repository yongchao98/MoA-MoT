def calculate_maximum_charge():
    """
    Calculates the maximum legal charge based on an estimate
    according to Ontario's Consumer Protection Act.
    """
    estimate = 3000.00
    percentage_increase_limit = 0.10  # 10%
    invoiced_amount = 3500.00

    # Calculate the maximum allowed increase and the total legal charge
    allowed_increase = estimate * percentage_increase_limit
    max_legal_charge = estimate + allowed_increase

    print("Legal Analysis based on Ontario's Consumer Protection Act:")
    print(f"Original Written Estimate: ${estimate:.2f}")
    print(f"Invoiced Amount: ${invoiced_amount:.2f}")
    print("\nThe Act states a supplier cannot charge more than 10% above the estimate.")
    print("\nCalculating the maximum legal charge:")
    # Printing the equation with all numbers as requested
    print(f"${estimate:.2f} (Estimate) + (${estimate:.2f} * {percentage_increase_limit:.2f}) = ${max_legal_charge:.2f}")
    print(f"\nThe maximum amount Marc is legally required to pay is ${max_legal_charge:.2f}.")
    print(f"The invoiced amount of ${invoiced_amount:.2f} is more than the legally permissible amount.")

calculate_maximum_charge()