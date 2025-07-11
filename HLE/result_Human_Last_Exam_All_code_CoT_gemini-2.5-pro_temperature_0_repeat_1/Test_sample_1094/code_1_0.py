import math

# This script calculates the normalized AC loss per cycle for a superconductor
# with an elliptical cross-section, based on the Norris model.

# The formula for the normalized loss 2*pi*Q/(mu_0*Ic^2) is:
# 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]
# where i = Im/Ic is the ratio of the current amplitude to the critical current.

print("--- AC Loss Calculator for an Elliptical Superconductor ---")
print("This script calculates the normalized loss: 2*pi*Q / (mu_0 * Ic^2)")
print("based on the normalized current i = Im/Ic.")
print("-" * 60)

try:
    # Prompt the user to enter the value of i
    i_str = input("Enter the value of the normalized current i (where 0 <= i < 1): ")
    i = float(i_str)

    # Validate that the input 'i' is within the allowed range [0, 1)
    if not (0 <= i < 1):
        print("\nError: The value of i must be greater than or equal to 0 and less than 1.")
        print("Please run the script again with a valid input.")
    else:
        # Handle the trivial case where i=0 (no current, no loss)
        if i == 0:
            term1 = 0.0
            term2 = 0.0
            term3 = 0.0
            normalized_loss = 0.0
        else:
            # Calculate each component of the Norris formula
            # term1 corresponds to (1-i)*ln(1-i)
            term1 = (1 - i) * math.log(1 - i)
            
            # term2 corresponds to (1+i)*ln(1+i)
            term2 = (1 + i) * math.log(1 + i)
            
            # term3 corresponds to i^2
            term3 = i**2
            
            # Combine the terms and multiply by the constant factor 2
            normalized_loss = 2 * (term1 + term2 - term3)

        # Print a detailed breakdown of the calculation
        print("\n--- Calculation Breakdown ---")
        print(f"The final equation is: Normalized Loss = 2 * ( term1 + term2 - term3 )")
        print(f"For your input i = {i}:")
        print(f"  - Constant Multiplier: 2")
        print(f"  - term1 = (1-i)*ln(1-i) = ({1-i:.4f})*ln({1-i:.4f}) = {term1:.6f}")
        print(f"  - term2 = (1+i)*ln(1+i) = ({1+i:.4f})*ln({1+i:.4f}) = {term2:.6f}")
        print(f"  - term3 = i^2 = {i:.4f}^2 = {term3:.6f}")
        
        print("\n--- Final Result ---")
        print(f"The normalized loss 2*pi*Q/(mu_0*Ic^2) is: {normalized_loss:.6f}")

except ValueError:
    print("\nError: Invalid input. Please enter a numerical value for i.")
except Exception as e:
    print(f"\nAn unexpected error occurred: {e}")
