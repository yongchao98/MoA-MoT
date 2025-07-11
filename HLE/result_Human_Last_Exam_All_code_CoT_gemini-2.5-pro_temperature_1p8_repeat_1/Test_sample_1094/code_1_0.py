import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for a superconducting elliptic bar.

    The normalized loss is given by 2*pi*Q / (mu_0 * Ic^2), where Q is the loss
    per cycle per unit length, Ic is the critical current, and i = Im/Ic is the
    normalized amplitude of the transport AC current.

    The calculation is based on the Norris formula, valid for i < 1.
    """
    # The formula is valid for 0 <= i < 1.
    if not (0 <= i < 1):
        print("Error: The normalized current 'i' must be in the range [0, 1).")
        return

    # Handle the edge case i=0, where loss is zero.
    if i == 0:
        result = 0.0
        print("For i = 0, the normalized loss is 0.")
        return

    # Calculate each term of the equation
    # L_norm = 2 * [ term1 + term2 - term3 ]
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2

    result = 2 * (term1 + term2 - term3)

    # Print the detailed breakdown as requested
    print(f"Calculating the normalized loss for i = {i}:")
    print("\nThe final equation for the normalized loss L_norm is:")
    print("L_norm = 2 * ((1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i^2)")

    print("\nPlugging in the numbers for each part of the equation:")
    print(f"  (1 - i)         = {1-i:.4f}")
    print(f"  ln(1 - i)       = {math.log(1-i):.4f}")
    print(f"  (1 - i)*ln(1-i) = {term1:.4f}")
    print("  --------------------")
    print(f"  (1 + i)         = {1+i:.4f}")
    print(f"  ln(1 + i)       = {math.log(1+i):.4f}")
    print(f"  (1 + i)*ln(1+i) = {term2:.4f}")
    print("  --------------------")
    print(f"  i^2             = {term3:.4f}")

    print("\nPutting it all together:")
    print(f"L_norm = 2 * ({term1:.4f} + {term2:.4f} - {term3:.4f})")
    print(f"L_norm = 2 * ({term1 + term2 - term3:.4f})")
    
    print(f"\nFinal Normalized Loss = {result}")

if __name__ == "__main__":
    try:
        # Prompt user for the normalized current 'i'
        i_input_str = input("Enter the value for the normalized current i = Im/Ic (e.g., 0.5): ")
        i_value = float(i_input_str)
        calculate_normalized_loss(i_value)
    except ValueError:
        print("Invalid input. Please enter a numerical value for i.")
    except Exception as e:
        print(f"An error occurred: {e}")