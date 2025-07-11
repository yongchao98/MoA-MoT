import math

def calculate_superconductor_ac_loss():
    """
    Calculates the normalized AC loss for a superconducting bar with an elliptical cross-section.
    """
    # Introduction and formula presentation
    print("This script calculates the normalized AC loss per cycle per unit length [2πQ/(μ₀*Ic²)]")
    print("for a superconducting elliptic bar carrying a transport AC current.")
    print("The calculation is valid for i = Im/Ic < 1, where Im is the current amplitude and Ic is the critical current.")
    print("\nThe formula derived by W.T. Norris is:")
    print("Loss = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i² ]\n")

    try:
        # Prompt user for the reduced current 'i'
        i_str = input("Enter the value for the reduced current i (a number between 0 and 1): ")
        i = float(i_str)

        # Validate the input 'i'
        if not (0 < i < 1):
            print("\nError: The value of i must be strictly between 0 and 1 for this formula to be valid.")
            return

        # --- Calculation steps ---
        # These intermediate variables are used to show the breakdown of the calculation.
        one_minus_i = 1 - i
        one_plus_i = 1 + i
        i_squared = i**2

        # Calculate the natural logarithm terms
        # math.log is the natural logarithm (ln)
        log_one_minus_i = math.log(one_minus_i)
        log_one_plus_i = math.log(one_plus_i)

        # Calculate the full terms inside the brackets
        term1 = one_minus_i * log_one_minus_i
        term2 = one_plus_i * log_one_plus_i

        # Calculate the final result
        sum_in_brackets = term1 + term2 - i_squared
        normalized_loss = 2 * sum_in_brackets

        # --- Output the detailed calculation ---
        print("\n--- Calculation Breakdown ---")
        print(f"For i = {i:.4f}:")
        
        # Display the equation with the user's numbers plugged in
        print("\n1. The equation with the numbers inserted is:")
        print(f"Loss = 2 * [ ({one_minus_i:.4f})*ln({one_minus_i:.4f}) + ({one_plus_i:.4f})*ln({one_plus_i:.4f}) - ({i:.4f})² ]")

        # Display the value of each term
        print("\n2. Evaluating each part of the equation:")
        print(f"Loss = 2 * [ ({term1:.6f}) + ({term2:.6f}) - ({i_squared:.6f}) ]")
        
        # Display the result of the expression in the brackets
        print("\n3. Simplifying the expression inside the brackets:")
        print(f"Loss = 2 * [ {sum_in_brackets:.6f} ]")
        
        # Display the final calculated loss
        print("\n4. Final Result:")
        print(f"The normalized loss 2πQ/(μ₀*Ic²) is: {normalized_loss:.6f}")

    except ValueError:
        print("\nError: Invalid input. Please enter a valid number.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

# Execute the function
calculate_superconductor_ac_loss()