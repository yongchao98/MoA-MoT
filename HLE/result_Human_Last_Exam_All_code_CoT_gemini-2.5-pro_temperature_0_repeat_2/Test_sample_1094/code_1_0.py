import math

def print_loss_formula():
    """
    This function prints the derived formula for the normalized AC loss
    in an elliptical superconductor carrying a transport current.

    The problem asks for the loss per cycle per unit length (Q) as a function
    of i = Im/Ic, presented in the standard form: 2*pi*Q/(mu_0*Ic^2).

    According to the critical state model (Norris formula), for i < 1, this
    normalized loss is independent of the ellipse's aspect ratio a/b.

    The final equation contains the numbers 2, 1, and the exponent 2, as requested.
    """

    # The derived formula for the normalized loss
    final_equation = "2 * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i**2 ]"

    print("The normalized loss per cycle, 2*pi*Q/(mu_0*Ic^2), as a function of i = Im/Ic is:")
    print(final_equation)

# Execute the function to print the result
print_loss_formula()

# For reference, a function to calculate the value for a given i
def calculate_normalized_loss(i):
    """
    Calculates the numerical value of the normalized AC loss.
    
    Args:
        i (float): The normalized current amplitude (Im/Ic). Must be between 0 and 1.

    Returns:
        float: The calculated normalized loss.
    """
    if not (0 <= i < 1):
        raise ValueError("The normalized current 'i' must be in the range [0, 1).")
    if i == 0:
        return 0.0
    
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    return 2 * (term1 + term2 - term3)
