import math

def calculate_normalized_ac_loss():
    """
    Calculates the normalized AC loss for a superconducting elliptic bar.

    This script calculates the loss per cycle per unit length (Q) for a
    superconducting elliptic bar in the critical state, carrying a transport
    AC current with amplitude Im, for the case where i = Im/Ic < 1.

    The result is given in the standard normalized form: 2*pi*Q / (mu_0 * Ic^2),
    which is calculated using the Norris formula:
    Loss(i) = 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2]
    """
    print("This script calculates the normalized AC loss for a superconducting elliptic bar.")
    print("The normalized loss is given by the function: 2*pi*Q / (mu_0 * Ic^2)")
    print("Formula: 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2]")
    print("-" * 50)

    while True:
        try:
            # Prompt the user for the normalized current i = Im/Ic
            i_str = input("Enter the value of i = Im/Ic (must be between 0 and 1, e.g., 0.5): ")
            i = float(i_str)

            # Validate the input range
            if not (0 <= i < 1):
                print("Error: The value of 'i' must be greater than or equal to 0 and strictly less than 1.")
                continue

            break # Exit loop if input is valid
        except ValueError:
            print("Error: Please enter a valid number.")

    # Handle the edge case i = 0 where ln(1-i) = ln(1) = 0
    if i == 0:
        normalized_loss = 0.0
        print("\nFor i = 0, the normalized loss is 0.")
    else:
        # Calculate each part of the equation
        val_1_minus_i = 1 - i
        val_1_plus_i = 1 + i
        val_i_squared = i**2

        # Use math.log for the natural logarithm (ln)
        ln_1_minus_i = math.log(val_1_minus_i)
        ln_1_plus_i = math.log(val_1_plus_i)
        
        term1 = val_1_minus_i * ln_1_minus_i
        term2 = val_1_plus_i * ln_1_plus_i
        term3 = val_i_squared

        # Calculate the final result
        sum_inside_brackets = term1 + term2 - term3
        normalized_loss = 2 * sum_inside_brackets

        # Print the detailed breakdown of the calculation
        print(f"\nCalculating for i = {i}:")
        print("-" * 25)
        print("Equation: 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2]")
        print("\n--- Substituting numbers ---")
        # Print each number in the equation
        print(f"1 - i = {val_1_minus_i}")
        print(f"1 + i = {val_1_plus_i}")
        print(f"i^2 = {val_i_squared:.6f}")
        print(f"ln(1-i) = ln({val_1_minus_i}) = {ln_1_minus_i:.6f}")
        print(f"ln(1+i) = ln({val_1_plus_i}) = {ln_1_plus_i:.6f}")
        print("\n--- Calculating terms ---")
        print(f"(1-i)*ln(1-i) = ({val_1_minus_i}) * ({ln_1_minus_i:.6f}) = {term1:.6f}")
        print(f"(1+i)*ln(1+i) = ({val_1_plus_i}) * ({ln_1_plus_i:.6f}) = {term2:.6f}")
        print("\n--- Final Calculation ---")
        print(f"Expression = 2 * [ {term1:.6f} + {term2:.6f} - {term3:.6f} ]")
        print(f"Expression = 2 * [ {sum_inside_brackets:.6f} ]")
        print(f"\nFinal Normalized Loss = {normalized_loss:.6f}")

if __name__ == '__main__':
    calculate_normalized_ac_loss()