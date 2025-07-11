import math

def calculate_normalized_elliptic_loss(i):
    """
    Calculates the normalized AC loss for a superconducting elliptic bar carrying
    a transport AC current, based on Norris's formula for the critical state model.

    The normalized loss is 2*pi*Q / (mu_0 * Ic^2), where Q is the loss per
    cycle per unit length, mu_0 is the permeability of free space, Ic is the
    critical current, and i = Im / Ic is the normalized current amplitude (i < 1).

    This normalized loss is notably independent of the ellipse's aspect ratio.

    Args:
        i (float): The normalized current amplitude, Im/Ic. Must be 0 < i < 1.
    """
    if not (0 < i < 1):
        print("Error: The normalized current 'i' must be greater than 0 and less than 1.")
        return

    # Norris's formula for normalized loss L = 2*pi*Q / (mu_0 * Ic^2) is:
    # L = 2 * [ (1-i) * ln(1-i) + (1+i) * ln(1+i) - i^2 ]
    
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    normalized_loss = 2 * (term1 + term2 - term3)

    print("The formula for the normalized AC loss 'L = 2*pi*Q / (mu_0*Ic^2)' as a function of 'i = Im/Ic' is:")
    print("L = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]\n")

    print(f"For a normalized current i = {i}, the calculation is as follows:")
    
    # Printing the equation with substituted numbers
    print("Step 1: Substitute i into the formula.")
    # Showing each number in the equation
    print(f"L = 2 * [ ({1-i})*ln({1-i}) + ({1+i})*ln({1+i}) - {i}^2 ]")
    
    print("\nStep 2: Calculate each term.")
    print(f"  Term (1-i)*ln(1-i) = ({1-i}) * ln({1-i}) = {term1:.6f}")
    print(f"  Term (1+i)*ln(1+i) = ({1+i}) * ln({1+i}) = {term2:.6f}")
    print(f"  Term i^2 = {i}^2 = {term3:.6f}")

    print("\nStep 3: Combine the terms and compute the final result.")
    print(f"L = 2 * [ {term1:.6f} + {term2:.6f} - {term3:.6f} ]")
    print(f"L = 2 * [ {term1 + term2 - term3:.6f} ]")
    print(f"L = {normalized_loss:.6f}\n")
    
    print("---------------------------------------------------------------")
    print(f"The final normalized loss for i = {i} is: {normalized_loss}")
    print("---------------------------------------------------------------")

# --- You can change this value for your own calculation ---
# The normalized current 'i' must be between 0 and 1 (e.g., 0.1, 0.5, 0.9).
normalized_current_i = 0.5
# ---------------------------------------------------------

# Execute the calculation
calculate_normalized_elliptic_loss(normalized_current_i)