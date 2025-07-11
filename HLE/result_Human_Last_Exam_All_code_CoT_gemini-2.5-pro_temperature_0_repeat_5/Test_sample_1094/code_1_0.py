import math

def calculate_normalized_loss(i):
    """
    Calculates and prints the normalized AC loss for an elliptical superconductor.

    The function calculates the loss for a given normalized current i = Im/Ic,
    where i must be between 0 and 1.

    The formula used is derived from Norris's model for elliptical wires:
    Normalized Loss = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]
    """
    if not (0 < i < 1):
        print("Error: The normalized current 'i' must be between 0 and 1.")
        return

    print(f"Calculating normalized AC loss for i = Im/Ic = {i}\n")

    # Print the general formula
    print("The general formula for the normalized loss 2*pi*Q / (mu_0 * Ic^2) is:")
    print("L_norm = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]\n")

    # Calculate each term
    val_1_minus_i = 1 - i
    val_1_plus_i = 1 + i
    val_i_squared = i**2

    log_1_minus_i = math.log(val_1_minus_i)
    log_1_plus_i = math.log(val_1_plus_i)

    term1 = val_1_minus_i * log_1_minus_i
    term2 = val_1_plus_i * log_1_plus_i
    term3 = val_i_squared

    # Print the breakdown of the calculation
    print("Step-by-step calculation:")
    print(f"1-i = {val_1_minus_i}")
    print(f"1+i = {val_1_plus_i}")
    print(f"i^2 = {val_i_squared:.4f}\n")

    print(f"ln(1-i) = ln({val_1_minus_i}) = {log_1_minus_i:.4f}")
    print(f"ln(1+i) = ln({val_1_plus_i}) = {log_1_plus_i:.4f}\n")

    print("The terms inside the bracket are:")
    print(f"(1-i)*ln(1-i) = {val_1_minus_i} * {log_1_minus_i:.4f} = {term1:.4f}")
    print(f"(1+i)*ln(1+i) = {val_1_plus_i} * {log_1_plus_i:.4f} = {term2:.4f}")
    print(f"i^2 = {term3:.4f}\n")

    # Final calculation
    bracket_term = term1 + term2 - term3
    normalized_loss = 2 * bracket_term

    print("Putting it all together in the equation:")
    print(f"L_norm = 2 * [ {term1:.4f} + {term2:.4f} - {term3:.4f} ]")
    print(f"L_norm = 2 * [ {bracket_term:.4f} ]")
    print(f"L_norm = {normalized_loss:.4f}\n")

    print(f"The final normalized loss for i = {i} is: {normalized_loss}")

# --- Main execution ---
# You can change this value to calculate the loss for a different normalized current.
# The value must be between 0 and 1.
example_i = 0.6

calculate_normalized_loss(example_i)