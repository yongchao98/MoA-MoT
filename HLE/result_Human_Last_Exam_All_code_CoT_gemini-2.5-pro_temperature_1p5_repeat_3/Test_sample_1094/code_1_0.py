import math

def calculate_normalized_ac_loss(i):
    """
    Calculates the normalized AC loss for an elliptical superconductor.

    The function computes 2*pi*Q/(mu_0*Ic^2) based on the Norris model
    for an elliptical wire carrying an AC transport current with i = Im/Ic < 1.

    Args:
        i (float): The normalized current amplitude (Im/Ic). Must be between 0 and 1.
    """
    # Check if the input value for i is valid for this formula.
    if not (0 < i < 1):
        print("Error: The value of i = Im/Ic must be strictly between 0 and 1 for this formula.")
        return

    # Define the components of the calculation
    one_minus_i = 1 - i
    one_plus_i = 1 + i
    i_squared = i**2
    
    # In Python's math library, math.log() is the natural logarithm (ln)
    log_one_minus_i = math.log(one_minus_i)
    log_one_plus_i = math.log(one_plus_i)
    
    term1_product = one_minus_i * log_one_minus_i
    term2_product = one_plus_i * log_one_plus_i
    
    # This is the expression inside the main brackets of the formula
    bracket_expression = term1_product + term2_product - i_squared
    
    # The final normalized loss is twice the bracket expression
    normalized_loss = 2 * bracket_expression
    
    # As requested, we will print the calculation step-by-step, showing each number in the equation.
    print(f"Calculating the normalized loss 2*pi*Q/(mu_0*Ic^2) for i = {i}")
    print("Formula: 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]")
    print("\n--- Calculation Steps ---")
    
    print("\n1. Substitute i into the equation's components:")
    print(f"   (1 - i) = (1 - {i}) = {one_minus_i}")
    print(f"   (1 + i) = (1 + {i}) = {one_plus_i}")
    print(f"   i^2 = {i}^2 = {i_squared}")

    print("\n2. Evaluate the natural logarithms (ln):")
    print(f"   ln({one_minus_i}) = {log_one_minus_i:.6f}")
    print(f"   ln({one_plus_i}) = {log_one_plus_i:.6f}")

    print("\n3. Calculate the full terms in the formula:")
    print(f"   (1-i)*ln(1-i) = {one_minus_i} * {log_one_minus_i:.6f} = {term1_product:.6f}")
    print(f"   (1+i)*ln(1+i) = {one_plus_i} * {log_one_plus_i:.6f} = {term2_product:.6f}")

    print("\n4. Evaluate the expression inside the brackets:")
    print(f"   [ {term1_product:.6f} + {term2_product:.6f} - {i_squared} ] = {bracket_expression:.6f}")

    print("\n5. Calculate the final normalized loss (multiply by 2):")
    print(f"   2 * {bracket_expression:.6f} = {normalized_loss:.6f}")

    # The final answer as a number for the specific i
    print(f"\nFor i = {i}, the normalized loss 2*pi*Q/(mu_0*Ic^2) is {normalized_loss:.6f}")


# --- Main Execution ---
# We will use an example value for the normalized current, i.
# For example, let's take a current amplitude that is 50% of the critical current.
i_value = 0.5
calculate_normalized_ac_loss(i_value)