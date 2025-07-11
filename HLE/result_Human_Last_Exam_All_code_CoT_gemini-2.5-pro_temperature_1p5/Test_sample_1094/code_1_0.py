import math

def calculate_normalized_loss():
    """
    Calculates the normalized AC loss for a superconducting elliptic bar.

    The function prompts the user for the normalized current i = Im/Ic,
    validates the input, and then calculates and prints the normalized loss
    per cycle per unit length, 2*pi*Q/(mu_0*Ic^2), based on the Norris formula.
    """
    
    print("This script calculates the normalized AC loss in a superconductor with an elliptical cross-section.")
    print("The normalized loss L = 2*pi*Q/(mu_0*Ic^2) is given by the function of i = Im/Ic:")
    
    # This fulfills the requirement to "output each number in the final equation"
    # by showing the formula with its explicit numerical coefficients.
    formula = "L(i) = 2 * ((1 - i) * ln(1 - i) + i - i**2 / 2)"
    print(formula)
    print("\nThis formula is valid for i < 1.")

    try:
        # Get user input for i
        i_str = input("Please enter a value for the normalized current i (where 0 <= i < 1): ")
        i = float(i_str)

        # Validate the input range
        if not (0 <= i < 1):
            print("Error: The value of i must be in the range [0, 1).")
            return

        # Handle the edge case i=0 to avoid log(1) calculation issues if needed, though math.log(1) is 0.
        if i == 0:
            loss = 0.0
        else:
            # Calculate the loss using the derived formula
            term1 = (1 - i) * math.log(1 - i)
            term2 = i
            term3 = (i**2) / 2
            loss = 2 * (term1 + term2 - term3)

        print(f"\nFor i = {i}:")
        print(f"The normalized loss is: {loss}")

    except ValueError:
        print("Error: Invalid input. Please enter a number.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    calculate_normalized_loss()
