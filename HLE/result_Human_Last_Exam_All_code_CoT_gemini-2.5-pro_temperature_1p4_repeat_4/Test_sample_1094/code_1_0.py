import math

def calculate_normalized_loss():
    """
    Calculates the normalized AC loss for an elliptic superconducting bar
    carrying a transport current with amplitude below the critical current.
    """
    # The normalized AC loss is a function of i = Im/Ic, where Im is the
    # current amplitude and Ic is the critical current.
    # The calculation is valid for i < 1.

    # We will use a representative value for i.
    # The user can modify this value to calculate the loss for a different current.
    i = 0.5

    print("This script calculates the normalized AC loss 2*pi*Q / (mu_0 * Ic^2) for an elliptic superconductor.")
    print("The analytical formula for the normalized loss as a function of i = Im/Ic is:")
    print("Loss(i) = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]\n")

    print(f"--- Calculation for i = {i} ---")

    # Check if i is in the valid range
    if not 0 < i < 1:
        print("Error: The value of i must be between 0 and 1 for this formula to be valid.")
        return

    # --- Step-by-step calculation ---

    # 1. Define all the components of the equation based on i
    val_1_minus_i = 1 - i
    val_1_plus_i = 1 + i
    val_i_squared = i**2

    # 2. Calculate the natural logarithms
    log_1_minus_i = math.log(val_1_minus_i)
    log_1_plus_i = math.log(val_1_plus_i)

    # 3. Calculate the individual terms
    term1 = val_1_minus_i * log_1_minus_i
    term2 = val_1_plus_i * log_1_plus_i

    # 4. Calculate the value inside the main brackets
    bracket_result = term1 + term2 - val_i_squared

    # 5. Calculate the final normalized loss
    normalized_loss = 2 * bracket_result

    # --- Output the results with each number in the equation ---
    print("Final Equation:")
    print(f"Loss = 2 * [ (1 - {i}) * ln(1 - {i}) + (1 + {i}) * ln(1 + {i}) - {i}^2 ]")
    print("\nSubstituting numerical values:")
    print(f"Loss = 2 * [ ({val_1_minus_i}) * ln({val_1_minus_i}) + ({val_1_plus_i}) * ln({val_1_plus_i}) - {val_i_squared} ]")
    print(f"Loss = 2 * [ ({val_1_minus_i}) * ({log_1_minus_i:.6f}) + ({val_1_plus_i}) * ({log_1_plus_i:.6f}) - {val_i_squared} ]")
    print(f"Loss = 2 * [ {term1:.6f} + {term2:.6f} - {val_i_squared} ]")
    print(f"Loss = 2 * [ {bracket_result:.6f} ]")
    print(f"Loss = {normalized_loss:.6f}\n")

    print("The final normalized loss per cycle is:")
    print(normalized_loss)


# Execute the function
if __name__ == "__main__":
    calculate_normalized_loss()
