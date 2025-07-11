import math

def calculate_maximum_chargeable_amount():
    """
    Calculates the maximum legal amount to be charged based on the Ontario Consumer Protection Act.
    """
    # Initial figures from the problem description
    estimated_price = 3000.00
    invoiced_price = 3500.00

    # Rule from the Ontario Consumer Protection Act
    allowable_overage_percentage = 0.10 # 10%

    # Calculation
    max_allowable_overage_amount = estimated_price * allowable_overage_percentage
    max_legal_price = estimated_price + max_allowable_overage_amount

    # Output the explanation and result
    print("Under the Ontario Consumer Protection Act, a business cannot charge more than 10% over a written estimate without the consumer's prior approval.")
    print("\nHere is the calculation based on the details provided:")
    print("---------------------------------------------------------")
    print(f"Original Written Estimate: ${estimated_price:.2f}")
    print(f"Allowable Overage Percentage: {allowable_overage_percentage * 100:.0f}%")
    print(f"Maximum Allowable Overage Amount: ${max_allowable_overage_amount:.2f}")
    print("\nFinal Calculation:")
    print(f"${estimated_price:.2f} (Estimate) + ${max_allowable_overage_amount:.2f} (10% Overage) = ${max_legal_price:.2f}")
    print("---------------------------------------------------------")
    print(f"\nThe invoice for ${invoiced_price:.2f} exceeds the maximum legal price.")
    print(f"Therefore, the lawyer would advise Marc that he is only required to pay ${max_legal_price:.2f}.")

# Execute the function to display the result
calculate_maximum_chargeable_amount()
<<<B>>>