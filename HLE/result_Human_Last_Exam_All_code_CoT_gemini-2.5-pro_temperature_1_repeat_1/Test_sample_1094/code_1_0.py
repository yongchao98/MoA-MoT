import math
import sys

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptic superconductor based on Norris's formula.
    The normalized loss is given by 2*pi*Q / (mu_0 * Ic^2).
    The formula is: 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2]
    
    Args:
        i (float): The normalized current amplitude (Im/Ic), where 0 < i < 1.
        
    Returns:
        float: The calculated normalized loss.
    """
    # The formula is only valid for i < 1. i > 0 is required for log.
    if not (0 < i < 1):
        raise ValueError("The normalized current 'i' must be strictly between 0 and 1.")

    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2

    # The full formula for the normalized loss
    result = 2 * (term1 + term2 - term3)
    return result

def main():
    """
    Main function to get user input and print the detailed calculation for the normalized loss.
    """
    print("This script calculates the normalized AC loss per cycle, 2*pi*Q/(mu_0*Ic^2),")
    print("for a superconducting bar with an elliptical cross-section.")
    print("The result depends only on the normalized current i = Im/Ic.")
    
    try:
        # Prompt the user for the value of i
        i_str = input("\nPlease enter a value for the normalized current i (e.g., 0.5): ")
        i_val = float(i_str)

        # Calculate the loss using the provided value
        loss = calculate_normalized_loss(i_val)

        # Unpack values for the detailed printout, fulfilling the request to
        # "output each number in the final equation"
        val_1_minus_i = 1 - i_val
        val_1_plus_i = 1 + i_val
        val_i_squared = i_val**2

        print("\n--- Calculation Breakdown ---")
        print(f"The general formula for the normalized loss is: 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]")
        print("\nPlugging in your value i =", i_val)
        print("The equation becomes:")
        # The following print statement shows the numbers within the equation structure
        print(f"2 * [ ({val_1_minus_i:.3f})*ln({val_1_minus_i:.3f}) + ({val_1_plus_i:.3f})*ln({val_1_plus_i:.3f}) - {i_val:.3f}^2 ]")
        print(f"= 2 * [ ({(1-i_val)*math.log(1-i_val):.4f}) + ({(1+i_val)*math.log(1+i_val):.4f}) - {val_i_squared:.4f} ]")
        print(f"= 2 * [ {((1-i_val)*math.log(1-i_val) + (1+i_val)*math.log(1+i_val) - val_i_squared):.4f} ]")
        
        print("\n-------------------------------------------------")
        print(f"The final normalized loss for i = {i_val} is: {loss}")
        print("-------------------------------------------------")

    except ValueError as ve:
        # Handle cases where input is not a number or out of the valid range 0 < i < 1
        print(f"\nError: {ve}. Please enter a numerical value strictly between 0 and 1.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()