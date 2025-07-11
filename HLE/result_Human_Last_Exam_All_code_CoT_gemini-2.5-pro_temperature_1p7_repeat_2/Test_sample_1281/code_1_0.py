def calculate_final_price():
    """
    Calculates the maximum legal price based on an original estimate and a
    percentage-based overage limit, as per Ontario's Consumer Protection Act.
    """
    original_estimate = 3000.00
    invoiced_price = 3500.00
    # The Consumer Protection Act allows for the final price to exceed the
    # estimate by no more than 10%.
    overage_percentage_limit = 0.10

    # Step 1: Calculate the maximum allowable dollar increase over the estimate.
    max_increase_amount = original_estimate * overage_percentage_limit

    # Step 2: Calculate the maximum total legal price.
    max_legal_price = original_estimate + max_increase_amount

    print("This script calculates the maximum legal amount Marc has to pay.")
    print("----------------------------------------------------------------")
    print(f"Original written estimate: ${original_estimate:.2f}")
    print(f"Final invoice amount: ${invoiced_price:.2f}")
    print(f"Maximum allowed increase percentage: {overage_percentage_limit * 100:.0f}%")
    print("\nCalculation Steps:")
    print("1. Determine the maximum allowed increase in dollars:")
    print(f"   ${original_estimate:.2f} (Estimate) * {overage_percentage_limit} (Limit) = ${max_increase_amount:.2f}")
    print("\n2. Determine the total maximum legal price:")
    print(f"   ${original_estimate:.2f} (Estimate) + ${max_increase_amount:.2f} (Max Increase) = ${max_legal_price:.2f}")
    print("----------------------------------------------------------------")
    print(f"\nConclusion: Marc is legally required to pay a maximum of ${max_legal_price:.2f}.")

calculate_final_price()
<<<B>>>