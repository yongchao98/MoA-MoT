import math

def calculate_norris_loss(i):
    """
    Calculates the normalized AC loss for an elliptical superconductor.

    The function computes the value of 2 * [ (1-i)ln(1-i) + (1+i)ln(1+i) - i^2 ],
    which represents the normalized loss 2*pi*Q / (mu_0 * Ic^2).

    Args:
        i (float): The normalized current amplitude (Im/Ic). Must be between 0 and 1.

    Returns:
        float: The calculated normalized loss, or None if the input is invalid.
    """
    # Validate the input 'i'
    if not (0 < i < 1):
        print(f"Error: The input value i={i} is invalid.")
        print("The formula is valid only for 0 < i < 1.")
        return None

    # --- Calculation Steps ---
    # The requirement is to show each number in the final equation.
    # We will calculate each part of the formula and print it.

    print("The general formula for the normalized loss is:")
    print("2 * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i^2 ]\n")
    
    print(f"For i = {i}, the equation becomes:")
    print(f"2 * [ (1 - {i})*ln(1 - {i}) + (1 + {i})*ln(1 + {i}) - {i}^2 ]\n")

    print("--- Step-by-Step Calculation ---")
    
    # Calculate each term
    term1_val = (1 - i) * math.log(1 - i)
    term2_val = (1 + i) * math.log(1 + i)
    term3_val = i**2
    
    print(f"Term 1: (1 - {i}) * ln(1 - {i}) = {1-i} * ln({1-i}) = {term1_val:.6f}")
    print(f"Term 2: (1 + {i}) * ln(1 + {i}) = {1+i} * ln({1+i}) = {term2_val:.6f}")
    print(f"Term 3: {i}^2 = {term3_val:.6f}")
    
    # Calculate the expression inside the brackets
    bracket_expression = term1_val + term2_val - term3_val
    print(f"\nExpression in brackets: {term1_val:.6f} + {term2_val:.6f} - {term3_val:.6f} = {bracket_expression:.6f}")

    # Final result
    final_loss = 2 * bracket_expression
    print(f"\nFinal Calculation: 2 * {bracket_expression:.6f} = {final_loss:.6f}")
    
    return final_loss

if __name__ == "__main__":
    # Set a sample value for the normalized current i = Im/Ic
    # This value must be > 0 and < 1.
    i_value = 0.5
    
    print(f"Calculating normalized loss for a transport current i = Im/Ic = {i_value}\n")
    
    # Run the calculation and print the final result
    normalized_loss = calculate_norris_loss(i_value)
    
    if normalized_loss is not None:
        print("\n----------------------------------------------------")
        print(f"The final normalized loss 2*pi*Q / (mu_0 * Ic^2) for i = {i_value} is: {normalized_loss:.6f}")
        print("----------------------------------------------------")