import math
import argparse

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for a superconductor with an elliptical
    cross-section carrying a transport current, based on the Norris model.

    The normalized loss is given by: 2*pi*Q / (mu_0 * Ic^2)

    Args:
        i (float): The ratio of the AC current amplitude to the critical current (Im/Ic).
                   Must be in the range 0 < i < 1.

    Returns:
        float: The value of the normalized loss. Returns None if i is out of range.
    """
    if not (0 < i < 1):
        print("Error: The current ratio 'i' must be a value between 0 and 1 (exclusive).")
        return None

    # The Norris formula for normalized loss
    # result = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]
    result = 2 * ((1 - i) * math.log(1 - i) + (1 + i) * math.log(1 + i) - i**2)
    return result

def main():
    """
    Main function to parse arguments and print the calculation.
    """
    parser = argparse.ArgumentParser(
        description='Calculate normalized AC loss in an elliptical superconductor.'
    )
    parser.add_argument(
        'i',
        type=float,
        help='Ratio of current amplitude to critical current, i = Im/Ic. Must be between 0 and 1.'
    )
    args = parser.parse_args()

    i_val = args.i
    loss = calculate_normalized_loss(i_val)

    if loss is not None:
        # Per the instructions, output the final equation with the numbers plugged in.
        # This shows "each number in the final equation".
        print("The formula for normalized loss is: 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]")
        print("\nFor i = {}:".format(i_val))

        # Create the string for the equation with numbers
        # We use {:.4f} to format the floating-point numbers for cleaner output
        term1_val_str = f"({1 - i_val:.4f})*ln({1 - i_val:.4f})"
        term2_val_str = f"({1 + i_val:.4f})*ln({1 + i_val:.4f})"
        term3_val_str = f"{i_val**2:.4f}"
        
        final_equation = f"Final Equation: 2 * [ {term1_val_str} + {term2_val_str} - {term3_val_str} ]"

        print(final_equation)
        print(f"Result = {loss}")

if __name__ == '__main__':
    main()
