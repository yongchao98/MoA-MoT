import sys

def calculate_consumer_protection_limit():
    """
    Calculates the maximum chargeable amount based on the Ontario Consumer Protection Act's 10% rule on estimates.
    """
    try:
        # Initial written estimate given by Honest Roofers Inc.
        estimate_price = 3000.00
        
        # The final invoiced price
        invoiced_price = 3500.00
        
        # The percentage by which the final price can exceed the estimate under the Consumer Protection Act
        allowable_increase_percentage = 0.10
        
        # Calculate the maximum allowable increase in dollars
        allowable_increase_amount = estimate_price * allowable_increase_percentage
        
        # Calculate the maximum legal price that can be charged
        max_legal_price = estimate_price + allowable_increase_amount
        
        print("Under Ontario's Consumer Protection Act, a final price cannot exceed the written estimate by more than 10%.")
        print("Let's calculate the maximum amount Marc is required to pay.")
        print("-" * 60)
        
        print(f"Original Written Estimate: ${estimate_price:,.2f}")
        print(f"Invoiced Amount: ${invoiced_price:,.2f}")
        print(f"Allowable Increase Percentage: {allowable_increase_percentage * 100}%")
        print(f"Maximum Allowable Increase Amount: ${allowable_increase_amount:,.2f}")
        
        print("\nFinal Calculation:")
        print(f"The maximum legal price is the original estimate plus the allowable increase:")
        print(f"${estimate_price:,.2f} (Estimate) + ${allowable_increase_amount:,.2f} (10% Increase) = ${max_legal_price:,.2f}")
        
        print(f"\nConclusion: Marc is legally required to pay HR no more than ${max_legal_price:,.2f}.")
        print("The invoiced amount of $3,500.00 is unenforceable.")
        
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

# Execute the function
calculate_consumer_protection_limit()
<<<B>>>