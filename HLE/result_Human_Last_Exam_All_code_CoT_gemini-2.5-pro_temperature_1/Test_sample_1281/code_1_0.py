def calculate_maximum_chargeable_amount():
    """
    Calculates the maximum amount a business can charge based on a written estimate
    under Ontario's Consumer Protection Act.
    """
    # Original written estimate provided by Honest Roofers Inc.
    estimate_price = 3000.00

    # The invoiced price after the job was completed.
    invoiced_price = 3500.00

    # Under Ontario's Consumer Protection Act, the final price cannot exceed
    # the estimate by more than 10%.
    allowed_increase_percentage = 0.10

    # Calculate the maximum allowed increase in dollars.
    max_increase_amount = estimate_price * allowed_increase_percentage

    # Calculate the total maximum legal price that can be charged.
    max_legal_price = estimate_price + max_increase_amount

    print(f"Original Estimated Price: ${estimate_price:.2f}")
    print(f"Invoiced Price: ${invoiced_price:.2f}")
    print(f"Allowable increase percentage under the Consumer Protection Act: {allowed_increase_percentage * 100}%")
    print("-" * 30)
    print("Calculating the maximum legally chargeable amount...")
    print(f"The maximum price is the estimate plus 10% of the estimate.")
    print(f"Final Equation: ${estimate_price:.2f} + (${estimate_price:.2f} * {allowed_increase_percentage}) = ${max_legal_price:.2f}")
    print("-" * 30)
    print(f"Therefore, the maximum amount Marc is required to pay is ${max_legal_price:.2f}.")

# Execute the function
calculate_maximum_chargeable_amount()