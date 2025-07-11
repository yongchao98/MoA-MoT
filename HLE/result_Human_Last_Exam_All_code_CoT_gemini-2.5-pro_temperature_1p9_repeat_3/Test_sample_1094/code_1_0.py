import math

def calculate_superconductor_loss():
    """
    Calculates the normalized AC loss for a superconducting elliptic bar.

    The function prompts the user for the normalized current i = Im/Ic,
    validates it, and then calculates the normalized loss per cycle per unit length,
    2*pi*Q/(mu_0*Ic^2), based on the Norris formula for i < 1.
    
    The final equation for the normalized loss is:
    L(i) = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]
    """
    try:
        # Prompt user for the normalized current i
        i_str = input("Enter the value of the normalized current i = Im/Ic (where 0 <= i < 1): ")
        i = float(i_str)

        # Validate the input range
        if not (0 <= i < 1):
            print("Error: The value of i must be between 0 (inclusive) and 1 (exclusive).")
            return

        # Handle the edge case i=0 where loss is zero
        if i == 0:
            print("\nFor i = 0, the loss is 0.")
            print("\nThe final calculated value for 2*pi*Q/(mu_0*Ic^2) is: 0.0")
            return

        # Calculate the components of the equation
        # Term1: (1 - i) * ln(1 - i)
        term1_factor1 = 1 - i
        term1_factor2 = math.log(1 - i) # ln in math
        term1_val = term1_factor1 * term1_factor2
        
        # Term2: (1 + i) * ln(1 + i)
        term2_factor1 = 1 + i
        term2_factor2 = math.log(1 + i)
        term2_val = term2_factor1 * term2_factor2
        
        # Term3: -i^2
        term3_val = i**2
        
        # Final result
        result = 2 * (term1_val + term2_val - term3_val)

        # Print the detailed calculation steps
        print(f"\nCalculating the normalized loss for i = {i}:")
        print("\nFormula: 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]")
        
        print(f"\nStep 1: Substitute i = {i} into the formula")
        print(f"2 * [ ({1-i:.4f})*ln({1-i:.4f}) + ({1+i:.4f})*ln({1+i:.4f}) - {i:.4f}^2 ]")
        
        print(f"\nStep 2: Evaluate each term in the brackets")
        print(f"[ {term1_val:.6f} + {term2_val:.6f} - {term3_val:.6f} ]")
        
        print(f"\nStep 3: Sum the terms inside the brackets")
        print(f"[ {term1_val + term2_val - term3_val:.6f} ]")
        
        print(f"\nStep 4: Multiply by 2 for the final result")
        print(f"2 * {term1_val + term2_val - term3_val:.6f} = {result:.6f}")

        print(f"\nThe final calculated value for 2*pi*Q/(mu_0*Ic^2) is: {result}")

    except ValueError:
        print("Error: Invalid input. Please enter a numerical value.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Execute the function
calculate_superconductor_loss()