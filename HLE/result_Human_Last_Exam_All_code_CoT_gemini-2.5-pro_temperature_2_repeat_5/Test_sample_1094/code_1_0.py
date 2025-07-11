import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptic superconductor based on Norris's formula.

    The normalized loss is given by 2 * pi * Q / (mu_0 * Ic^2),
    where Q is the loss per cycle per unit length, and i = Im / Ic.

    Args:
        i (float): The ratio of the current amplitude to the critical current (Im/Ic).
                   Must be in the range 0 <= i < 1.

    Returns:
        float: The dimensionless normalized loss. Returns 0 for i=0,
               or raises a ValueError for i >= 1.
    """
    if not 0 <= i < 1:
        raise ValueError(f"Input i={i} is outside the valid range [0, 1).")

    # At i=0, the loss is 0. (1-0)*ln(1-0) + (1+0)*ln(1+0) - 0 = 0.
    if i == 0:
        return 0.0
    
    # Norris's formula for an elliptic bar (independent of aspect ratio).
    # The final expression for normalized loss F(i) = 2*pi*Q / (mu_0 * Ic^2) is:
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    normalized_loss = 2 * (term1 + term2 - term3)
    
    return normalized_loss

if __name__ == "__main__":
    print("This script calculates the normalized AC loss in a superconducting elliptic bar.")
    print("The calculation is based on Norris's analytical solution for i = Im/Ic < 1.")
    
    # Per the instructions, printing the final equation with all its numbers.
    print("\nThe final analytical function for the normalized loss F(i) = 2*pi*Q / (mu_0 * Ic^2) is:")
    print("F(i) = 2 * ( (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i**2 )")
    
    print("\n--- Example Calculations ---")
    print(f"{'i = Im/Ic':<12} | {'Normalized Loss':<20}")
    print("-" * 35)
    
    # Calculate and print for a few representative values of i
    try:
        for i_val in [0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
            loss = calculate_normalized_loss(i_val)
            print(f"{i_val:<12.2f} | {loss:<20.6f}")
    except ValueError as e:
        print(e)