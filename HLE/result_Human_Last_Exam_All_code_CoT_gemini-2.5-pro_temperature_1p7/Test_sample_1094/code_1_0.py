import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptical superconductor.

    The normalized loss is given by 2 * pi * Q / (mu_0 * Ic^2), where:
    - Q is the loss per cycle per unit length.
    - i = Im / Ic is the normalized current amplitude (Im < Ic).
    - mu_0 is the permeability of free space.
    - Ic is the critical current.

    Args:
        i (float): The normalized current amplitude, must be in the range [0, 1).

    Returns:
        float: The calculated normalized loss. Returns -1.0 on invalid input.
    """
    if not (0 <= i < 1):
        print(f"Error: Input 'i' must be in the range [0, 1). Received i = {i}")
        return -1.0
    
    # Handle the edge case i=0 to avoid math domain errors if floating point issues arise
    if i == 0:
        return 0.0

    # The formula is 2 * [(1+i)ln(1+i) + (1-i)ln(1-i) - i^2]
    term1 = (1 + i) * math.log(1 + i)
    term2 = (1 - i) * math.log(1 - i)
    term3 = i**2
    
    normalized_loss = 2 * (term1 + term2 - term3)
    
    return normalized_loss

# Main part of the script
if __name__ == "__main__":
    # First, print the derived equation as requested.
    # The numbers in the equation are 2, 1, 1, 1, 2 (from i^2)
    print("The formula for the normalized loss 2*pi*Q/(mu_0*Ic^2) as a function of i = Im/Ic is:")
    print("Normalized Loss = 2 * ((1+i)*ln(1+i) + (1-i)*ln(1-i) - i**2)")
    print("-" * 50)
    print("Calculating normalized loss for some example values of i:")

    # Example values for i
    example_i_values = [0.1, 0.3, 0.5, 0.7, 0.9, 0.99]

    for i_val in example_i_values:
        loss = calculate_normalized_loss(i_val)
        if loss != -1.0:
            print(f"For i = {i_val:.2f}, the normalized loss is: {loss:.6f}")
