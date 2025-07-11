def calculate_final_price():
    """
    Calculates the maximum legal price based on an estimate according to
    Ontario's Consumer Protection Act.
    """
    estimate = 3000.00
    invoiced_price = 3500.00
    allowable_increase_percentage = 0.10

    # Calculate the maximum allowable increase in dollars
    max_increase_amount = estimate * allowable_increase_percentage

    # Calculate the maximum legal price
    max_legal_price = estimate + max_increase_amount

    print("This script calculates the maximum amount Marc is required to pay under Ontario's Consumer Protection Act.")
    print("-" * 80)
    print(f"Original Written Estimate: ${estimate:.2f}")
    print(f"Allowable Increase Percentage: {allowable_increase_percentage * 100:.0f}%")
    print("\nCalculating the maximum legal price:")
    print(f"Maximum Increase Amount = ${estimate:.2f} (Estimate) * {allowable_increase_percentage} = ${max_increase_amount:.2f}")
    print(f"Maximum Legal Price = ${estimate:.2f} (Estimate) + ${max_increase_amount:.2f} (Increase) = ${max_legal_price:.2f}")
    print("-" * 80)
    print(f"The invoiced price was ${invoiced_price:.2f}, which is above the legal limit.")
    print(f"Therefore, Marc is only required to pay HR ${max_legal_price:.2f}.")

calculate_final_price()