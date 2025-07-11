import math

def calculate_norris_loss(i):
    """
    Calculates the normalized AC transport current loss for a superconductor
    with an elliptical cross-section based on the Norris model.

    Args:
        i (float): The ratio of the current amplitude to the critical current (Im/Ic).
                   Must be in the range [0, 1).

    Returns:
        float: The normalized loss per cycle, 2*pi*Q / (mu_0 * Ic^2).
    """
    if not 0 <= i < 1:
        raise ValueError("The normalized current 'i' must be in the range [0, 1).")
    
    # Handle the edge case i=0 to avoid log(1) computation if not necessary,
    # though math.log(1) is 0.
    if i == 0:
        return 0.0

    # Norris formula for normalized loss
    # f(i) = (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2
    # Normalized loss = 2 * f(i)
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    normalized_loss = 2 * (term1 + term2 - term3)
    
    return normalized_loss

if __name__ == "__main__":
    # The normalized loss is a function of i = Im/Ic.
    # The final equation for the normalized loss 2*pi*Q / (mu_0 * Ic^2) is:
    # 2 * ((1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i^2)
    
    # We can represent the numbers and operators in the equation as follows:
    factor = 2
    op1 = "(1 - i) * ln(1 - i)"
    op2 = "(1 + i) * ln(1 + i)"
    op3 = "i^2" # where the exponent is 2

    print("The formula for the normalized AC loss, 2*pi*Q / (mu_0 * Ic^2), as a function of i = Im/Ic is:")
    print(f"{factor} * ({op1} + {op2} - {op3})\n")
    
    # Let's calculate the value for a specific example, i = 0.5
    example_i = 0.5
    try:
        loss_value = calculate_norris_loss(example_i)
        print(f"For an example transport current with amplitude i = Im/Ic = {example_i}:")
        print(f"The normalized loss is: {loss_value}")
        
        # This is the final numerical answer for the example case i = 0.5
        # The format <<<value>>> is for automated checking.
        print(f"\n<<<{loss_value}>>>")

    except ValueError as e:
        print(f"Error: {e}")