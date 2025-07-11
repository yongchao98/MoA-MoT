import sys
import io

# A text-based chatbot, written in Python, that can inform the user about the laws and regulations in their jurisdiction.

def calculate_max_charge(estimate_price, percentage_limit):
    """
    Calculates the maximum legal charge based on an estimate and a percentage limit.

    Args:
        estimate_price (float): The initial estimated price.
        percentage_limit (float): The decimal representation of the percentage limit (e.g., 0.10 for 10%).

    Returns:
        float: The maximum legally allowed price.
    """
    # Calculate the maximum allowed increase in dollar amount
    allowed_increase = estimate_price * percentage_limit
    
    # Calculate the total maximum legal price
    max_legal_price = estimate_price + allowed_increase
    
    return max_legal_price, allowed_increase

def main():
    """
    Main function to solve the user's problem.
    """
    # Given values from the problem
    estimate_price = 3000.00
    final_invoice = 3500.00
    # The percentage limit set by the Ontario Consumer Protection Act
    percentage_limit = 0.10

    # Calculate the maximum legal price Marc has to pay
    max_legal_price, allowed_increase = calculate_max_charge(estimate_price, percentage_limit)

    # Print the explanation and the result
    print("Under Ontario's Consumer Protection Act, when a written estimate is provided, the final price cannot be more than 10% higher.")
    print("\nHere is the calculation based on the rules:")
    print(f"1. Original Estimate Provided: ${estimate_price:.2f}")
    print(f"2. Maximum Allowed Increase: ${estimate_price:.2f} * {percentage_limit:.2f} = ${allowed_increase:.2f}")
    print("\nTherefore, the maximum legal amount Marc can be charged is:")
    # This line outputs each number in the final equation
    print(f"${estimate_price:.2f} (Estimate) + ${allowed_increase:.2f} (10% Overage) = ${max_legal_price:.2f}")
    
    print(f"\nSince the final invoice of ${final_invoice:.2f} is higher than this amount, Marc is only legally required to pay ${max_legal_price:.2f}.")

if __name__ == '__main__':
    main()