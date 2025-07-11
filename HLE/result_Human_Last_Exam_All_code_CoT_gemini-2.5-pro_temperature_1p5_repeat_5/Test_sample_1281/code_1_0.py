import sys

def calculate_consumer_rights():
    """
    Calculates the maximum legal price based on Ontario's Consumer Protection Act
    and explains the reasoning for Marc's situation.
    """
    estimate = 3000.00
    invoiced_amount = 3500.00
    allowed_overage_percent = 0.10

    # Calculate the maximum allowable increase and total price
    max_increase_amount = estimate * allowed_overage_percent
    max_legal_price = estimate + max_increase_amount

    # Print the explanation
    print("According to Ontario's Consumer Protection Act, a business cannot charge a consumer more than 10% above a given estimate.")
    print("This law applies as the business (Honest Roofers Inc.) is located in Ontario, where the agreement was made.")
    print("\n--- Calculation ---")
    print(f"Original Estimate: ${estimate:,.2f}")
    # Here we show each number in the final equation
    print(f"Maximum Allowed Increase (10%): ${max_increase_amount:,.2f}")
    print(f"Final Equation for Maximum Price: ${estimate:,.2f} + ${max_increase_amount:,.2f} = ${max_legal_price:,.2f}")
    print("-------------------")
    print(f"\nThe invoiced price was ${invoiced_amount:,.2f}.")
    print(f"Since this is more than the maximum legally allowed price of ${max_legal_price:,.2f}, Marc is only legally required to pay ${max_legal_price:,.2f}.")
    print("This aligns with Answer B.")

if __name__ == "__main__":
    calculate_consumer_rights()
<<<B>>>