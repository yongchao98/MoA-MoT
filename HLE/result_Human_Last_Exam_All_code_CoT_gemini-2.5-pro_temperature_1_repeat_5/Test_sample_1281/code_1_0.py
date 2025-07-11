def calculate_final_bill():
    """
    Calculates the maximum payable amount based on an original estimate
    and a legally mandated maximum percentage increase.
    """
    original_estimate = 3000.00
    invoiced_amount = 3500.00
    max_increase_percent = 10.0

    # Calculate the maximum amount that can be charged
    increase_multiplier = 1 + (max_increase_percent / 100.0)
    max_chargeable_amount = original_estimate * increase_multiplier

    # The amount Marc is required to pay is the lesser of the invoiced amount
    # and the maximum legally allowed amount.
    final_payable_amount = min(invoiced_amount, max_chargeable_amount)

    print("Under the Ontario Consumer Protection Act, a final price cannot exceed a written estimate by more than 10%.")
    print("Let's calculate the maximum amount Marc is required to pay.")
    print("-" * 50)
    print(f"Original Estimate: ${original_estimate:.2f}")
    print(f"Maximum Allowed Increase: {int(max_increase_percent)}%")
    print(f"Calculation of Maximum Chargeable Amount:")
    print(f"${original_estimate:.2f} * (1 + {int(max_increase_percent)} / 100) = ${max_chargeable_amount:.2f}")
    print("-" * 50)
    print(f"Since the invoiced amount of ${invoiced_amount:.2f} exceeds the maximum allowed, Marc is required to pay ${final_payable_amount:.2f}.")

calculate_final_bill()