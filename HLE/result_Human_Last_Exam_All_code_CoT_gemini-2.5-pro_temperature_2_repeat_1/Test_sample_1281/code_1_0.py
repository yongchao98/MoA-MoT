def calculate_roofer_payment():
    """
    Calculates the legally required payment based on Ontario's Consumer Protection Act.

    The Act stipulates that the final price cannot exceed a written estimate by more than 10%.
    """
    # Define the values from the scenario
    estimate_price = 3000.00
    invoiced_price = 3500.00
    allowed_overage_rate = 0.10

    # Calculate the maximum allowed increase and the total maximum payable price
    allowed_increase_amount = estimate_price * allowed_overage_rate
    max_payable_price = estimate_price + allowed_increase_amount

    # Print the explanation and the step-by-step calculation
    print("The lawyer's advice would be based on Ontario's Consumer Protection Act.")
    print("The Act limits the final charge to 10% more than the written estimate.")
    print("-" * 50)
    print(f"Original Written Estimate: ${estimate_price:.2f}")
    print(f"Invoiced Amount: ${invoiced_price:.2f}")
    print("-" * 50)
    print("Calculating the maximum legal amount Marc has to pay:")
    
    # Per the instructions, printing each number in the final equation
    print(f"Equation: ${estimate_price:.2f} (Estimate) + (${estimate_price:.2f} * {allowed_overage_rate}) = ${max_payable_price:.2f}")
    print(f"This simplifies to: ${estimate_price:.2f} (Estimate) + ${allowed_increase_amount:.2f} (10% Increase) = ${max_payable_price:.2f}")
    print("-" * 50)

    # State the final conclusion
    print(f"The invoiced price of ${invoiced_price:.2f} is more than the maximum allowed price of ${max_payable_price:.2f}.")
    print(f"Therefore, Marc is only legally required to pay HR ${max_payable_price:.2f}.")

# Run the calculation and print the results
calculate_roofer_payment()
<<<B>>>