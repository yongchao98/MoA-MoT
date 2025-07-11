import math

def calculate_normalized_superconductor_loss(i):
    """
    Calculates and prints the normalized AC loss for a superconducting
    elliptic bar carrying a transport current, for i < 1.

    Args:
        i (float): The normalized current amplitude (Im/Ic), where Im is the
                   current amplitude and Ic is the critical current.
                   The value of i must be between 0 and 1.
    """
    if not (0 < i < 1):
        print("Error: The normalized current 'i' must be in the range 0 < i < 1.")
        return

    # The problem asks for the normalized loss per cycle per unit length,
    # given by 2*pi*Q/(mu_0*Ic^2).
    # This is calculated using the Norris formula for an elliptical conductor.
    # Formula: F(i) = (1+i)*ln(1+i) + (1-i)*ln(1-i) - i^2
    
    val_1_plus_i = 1 + i
    val_1_minus_i = 1 - i
    val_i_squared = i**2

    term1 = val_1_plus_i * math.log(val_1_plus_i)
    term2 = val_1_minus_i * math.log(val_1_minus_i)
    normalized_loss = term1 + term2 - val_i_squared

    # --- Output ---

    # 1. Print the general formula
    print("The normalized loss per cycle per unit length, 2*pi*Q/(mu_0*Ic^2), as a function of i = Im/Ic is given by the formula:")
    print("F(i) = (1+i)*ln(1+i) + (1-i)*ln(1-i) - i^2\n")

    # 2. Print the formula with the specific numbers for the given i
    print(f"For the specific case where i = {i}:")
    print("The equation is evaluated as follows:")
    
    # This step fulfills the requirement: "output each number in the final equation!"
    print(f"F({i}) = (1+{i})*ln(1+{i}) + (1-{i})*ln(1-{i}) - {i}^2")
    
    # 3. Print the final result
    print(f"\nThe calculated value for the normalized loss is: {normalized_loss}")


# --- Main Execution ---
# We will demonstrate the calculation for a sample value of i = 0.5.
# You can change this value to calculate the loss for any other i < 1.
example_i = 0.5
calculate_normalized_superconductor_loss(example_i)

<<<(1+i)*ln(1+i) + (1-i)*ln(1-i) - i**2>>>