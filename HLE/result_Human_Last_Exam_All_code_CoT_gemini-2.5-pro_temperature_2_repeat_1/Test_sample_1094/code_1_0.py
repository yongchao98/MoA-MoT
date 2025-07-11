import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptic superconductor.

    This function implements the Norris formula for the transport-current AC loss
    in a solid superconductor with an elliptic cross-section for i < 1.

    The result is given in the standard normalized form: 2 * pi * Q / (mu_0 * Ic^2),
    where Q is the loss per cycle per unit length, and i = Im / Ic.

    Args:
        i (float): The normalized current amplitude (Im/Ic). Must be between 0 and 1.

    Returns:
        float: The calculated normalized loss. Returns None if i is out of range.
    """
    if not (0 <= i < 1):
        print("Error: The value of i must be in the range [0, 1).")
        return None

    # Handle the edge case i=0 where ln(1) = 0 to avoid precision issues.
    if i == 0:
        return 0.0

    # The loss function F(i) = (1-i)ln(1-i) + (1+i)ln(1+i) - i^2
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    f_i = term1 + term2 - i**2
    
    # The normalized loss is 2 * F(i)
    normalized_loss = 2 * f_i
    
    return normalized_loss

# The final formula as a function of i = Im/Ic
# The equation is: 2 * ((1-i)*ln(1-i) + (1+i)*ln(1+i) - i**2)
print("The normalized AC loss, 2*pi*Q / (mu_0 * Ic^2), is given by the function:")
# The following print statement displays the final formula with all its numbers and parts.
print("f(i) = 2 * ((1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i**2)")
print("-" * 30)

# --- Example Calculation ---
# Let's calculate the normalized loss for a sample value, e.g., i = 0.5
example_i = 0.5
loss_value = calculate_normalized_loss(example_i)

if loss_value is not None:
    print(f"For i = Im/Ic = {example_i}, the normalized loss is: {loss_value:.6f}")
