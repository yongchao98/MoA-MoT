import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss 2*pi*Q / (mu_0 * Ic^2) for an 
    elliptical superconductor carrying a transport current.

    This formula, derived by W.T. Norris, is valid for the sub-critical 
    regime, i.e., i = Im/Ic < 1.

    Args:
        i (float): The ratio of current amplitude to critical current, Im/Ic.

    Returns:
        float: The value of the normalized loss.
    """
    # The formula is not valid for i >= 1 (full penetration).
    # It is also not defined for i > 1 due to ln(1-i).
    if not (0 <= i < 1):
        raise ValueError("The normalized current 'i' must be in the range [0, 1).")
    
    # Handle the edge case i=0, where loss is zero.
    if i == 0:
        return 0.0

    # Calculate the terms of the formula
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    # The complete formula for the normalized loss
    result = 2 * (term1 + term2 - term3)
    return result

# --- Main execution ---

# The formula for the normalized loss is presented below.
formula_str = "2 * ((1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i**2)"

print("This script provides the formula for the normalized AC loss in an elliptical superconductor.")
print("The normalized loss is given by 2*pi*Q / (mu_0 * Ic^2), as a function of i = Im/Ic (for i < 1).")
print("\nThe complete formula is:")
print(f"  f(i) = {formula_str}\n")

print("Here are the individual numbers and components of the final equation:")
print("  - The constant multiplicative factor is: 2")
print("  - The first term inside the parentheses is: (1 - i) * ln(1 - i)")
print("  - The second term inside the parentheses is: (1 + i) * ln(1 + i)")
print("  - The third term inside the parentheses is: -i**2\n")

# You can use the function to calculate the loss for a specific 'i'.
# Here is an example for i = 0.5
try:
    i_example = 0.5
    normalized_loss = calculate_normalized_loss(i_example)
    print(f"Example calculation for i = {i_example}:")
    print(f"The normalized loss is: {normalized_loss}")
except ValueError as e:
    print(f"Error during calculation: {e}")
