import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for a superconductor in the critical state
    carrying a transport current.

    The formula is based on the Norris model for an isolated conductor with a
    convex cross-section (e.g., ellipse, circle).

    Args:
        i (float): The ratio of the current amplitude to the critical current (Im/Ic).
                   Must be in the range 0 <= i < 1.

    Returns:
        float: The normalized loss per cycle, 2*pi*Q / (mu_0 * Ic^2).
    """
    if not 0 <= i < 1:
        raise ValueError("The normalized current 'i' must be in the range [0, 1).")

    # The formula returns 0 for i=0.
    if i == 0:
        return 0.0

    # The Norris formula for normalized loss: 2 * [(1-i)ln(1-i) + (1+i)ln(1+i) - i^2]
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    result = 2 * (term1 + term2 - term3)
    return result

def main():
    """
    Main function to demonstrate the AC loss calculation for various currents.
    """
    print("This script calculates the normalized AC loss in a superconducting elliptic bar.")
    print("The normalized loss is 2*pi*Q / (mu_0 * Ic^2), where i = Im / Ic.")
    print("The formula used is: 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2]\n")
    
    # Example values for i = Im/Ic
    example_currents = [0.1, 0.25, 0.5, 0.75, 0.9, 0.95]

    for i in example_currents:
        try:
            # Calculate the loss
            loss = calculate_normalized_loss(i)
            
            # As requested, output the equation with the numbers plugged in.
            print(f"For i = {i}:")
            print(f"2 * ( (1-{i})*ln(1-{i}) + (1+{i})*ln(1+{i}) - {i}**2 ) = {loss:.6f}")
            print("-" * 50)

        except ValueError as e:
            print(f"Error calculating for i = {i}: {e}")

if __name__ == "__main__":
    main()