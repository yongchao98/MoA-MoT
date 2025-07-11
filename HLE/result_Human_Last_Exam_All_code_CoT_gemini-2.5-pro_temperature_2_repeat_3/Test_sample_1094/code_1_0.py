import math

def calculate_normalized_ac_loss(i):
    """
    Calculates the normalized AC loss for a superconducting elliptic bar based on the Norris model.
    The normalized loss is given by 2 * pi * Q / (mu_0 * Ic^2).

    Args:
        i (float): The ratio of the current amplitude to the critical current (Im/Ic).
                   Must be in the range [0, 1).
    """
    # 1. Input Validation
    if not (0 <= i < 1):
        print(f"Error: The value of i = {i} is invalid. It must be in the range 0 <= i < 1.")
        return None

    # 2. State the formula being used
    # The normalized loss is derived from the Norris formula for an elliptical wire:
    # Q = (μ_0 * Ic^2 / π) * F(i)
    # where F(i) = (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2
    # The required normalization is 2*π*Q/(μ_0*Ic^2) = 2*π/(μ_0*Ic^2) * [(μ_0 * Ic^2 / π) * F(i)] = 2*F(i)
    
    print("This script calculates the normalized AC loss 2*π*Q/(μ_0*Ic^2) for a superconducting elliptic bar.")
    print(f"The calculation is performed for i = Im/Ic = {i}\n")

    print("The final formula for the normalized loss is:")
    print("Normalized Loss = 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2]\n")

    # 3. Perform the calculation and show the steps
    print("--- Calculation Steps ---")
    
    # Define terms for clarity and to show the numbers
    one_minus_i = 1 - i
    one_plus_i = 1 + i
    i_squared = i**2

    # Handle the edge case i=0 to avoid log(1) calculations for clarity, though mathematically fine.
    if i == 0.0:
        term1 = 0.0
    else:
        term1 = one_minus_i * math.log(one_minus_i)
    
    term2 = one_plus_i * math.log(one_plus_i)

    # Full expression in brackets
    f_i = term1 + term2 - i_squared

    # Final normalized loss
    normalized_loss = 2 * f_i
    
    # Print the equation with the specific numbers plugged in
    print("Plugging in the numbers into the formula:")
    # Rounding the intermediate values for a cleaner display of the equation itself
    equation_str = (
        f"Normalized Loss = 2 * [({1-i:.2f})*ln({1-i:.2f}) + "
        f"({1+i:.2f})*ln({1+i:.2f}) - {i:.2f}^2]"
    )
    print(equation_str)
    
    # Show the value of each term
    print(f"Calculated terms: (1-i)*ln(1-i) = {term1}, (1+i)*ln(1+i) = {term2}, i^2 = {i_squared}")
    
    # Show the result of the expression in brackets
    print(f"Normalized Loss = 2 * [{f_i}]")
    
    # 4. Print the final result
    print("\n-------------------------------------------")
    print(f"The final normalized loss is: {normalized_loss}")
    print("-------------------------------------------")

    return normalized_loss

if __name__ == "__main__":
    # A representative value for i = Im/Ic, where i < 1
    i_value = 0.5
    
    # Execute the calculation
    final_answer = calculate_normalized_ac_loss(i_value)
    # The final answer for the <<<...>>> format is printed below,
    # it corresponds to the calculated loss for i=0.5.
    if final_answer is not None:
        # Formatted to prevent scientific notation in the final output block.
        final_answer_string = f"<<<{final_answer:.17f}>>>"
        # The line below is for the final answer block and will be stripped out,
        # but let's calculate and hold it.
        # print(final_answer_string) 
