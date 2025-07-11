import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptical superconductor in the low-current regime.
    The loss is given by the Norris formula, normalized as 2*pi*Q / (mu_0 * Ic^2).
    
    Args:
        i (float): The ratio of current amplitude to critical current, i = Im/Ic.
                   The value must be in the range 0 <= i < 1.

    Returns:
        float: The normalized loss, a dimensionless quantity.
    """
    if not 0 <= i < 1:
        # The formula is only valid for i < 1.
        raise ValueError("Input 'i' must be in the range [0, 1).")

    # Handle the edge case i=0, where loss is zero.
    if i == 0.0:
        return 0.0
    
    # The core of the Norris formula for i < 1
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    # The final normalized expression is 2 * ( ... )
    normalized_loss = 2 * (term1 + term2 - term3)
    
    return normalized_loss

if __name__ == "__main__":
    # First, print the general formula that was derived.
    print("This script calculates the AC loss in an elliptical superconductor for i = Im/Ic < 1.")
    print("The result is provided in the standard normalized form:")
    # Using 'ln' for natural logarithm to be mathematically precise.
    print("\n  2*pi*Q / (mu_0 * Ic^2) = 2 * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i^2 ]\n")
    print("-" * 70)
    print("Below are example calculations using this formula.")
    print("-" * 70)
    
    # Demonstrate the calculation for a few values of i.
    example_current_ratios = [0.25, 0.5, 0.75, 0.95]
    
    for i_val in example_current_ratios:
        # Calculate the result using the function
        loss_result = calculate_normalized_loss(i_val)
        
        print(f"For a current ratio i = {i_val}:")
        
        # Output the equation with the specific numbers substituted, as requested.
        # The numbers in the equation are the constant '2' and the variable 'i'.
        print(f"  Equation with numbers: Result = 2 * [ (1 - {i_val})*ln(1 - {i_val}) + (1 + {i_val})*ln(1 + {i_val}) - {i_val}^2 ]")
        
        # Print the final calculated value
        print(f"  Calculated Normalized Loss = {loss_result:.8f}\n")
