def calculate_roofing_cost():
    """
    Calculates the maximum legal price based on Ontario's Consumer Protection Act.
    
    The Act stipulates that a final price cannot exceed a written estimate by more than 10%.
    """
    
    # The initial written estimate given to Marc.
    estimate_price = 3000.00
    
    # The percentage by which the estimate can be legally exceeded.
    allowed_overage_percentage = 0.10 # This represents 10%
    
    # The invoiced amount from Honest Roofers Inc.
    invoiced_price = 3500.00
    
    # Calculate the dollar value of the allowed increase.
    allowed_increase_amount = estimate_price * allowed_overage_percentage
    
    # Calculate the maximum price that can be legally charged.
    max_legal_price = estimate_price + allowed_increase_amount
    
    print("Analyzing the transaction under Ontario's Consumer Protection Act:")
    print("-" * 60)
    print(f"Original written estimate: ${estimate_price:.2f}")
    print(f"Invoiced amount: ${invoiced_price:.2f}")
    print(f"Maximum allowed increase over estimate: {int(allowed_overage_percentage * 100)}%")
    print("-" * 60)
    print("The final amount Marc is legally required to pay is calculated as follows:")
    
    # Outputting the final equation with each number as requested.
    print(f"Maximum Legal Price = ${estimate_price:.2f} + ${allowed_increase_amount:.2f} = ${max_legal_price:.2f}")

calculate_roofing_cost()