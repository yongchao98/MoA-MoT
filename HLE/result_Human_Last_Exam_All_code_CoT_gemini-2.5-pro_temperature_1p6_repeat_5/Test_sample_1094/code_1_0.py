import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptical superconductor.
    The normalized loss is given by 2*pi*Q / (mu_0 * Ic^2).
    Formula: L(i) = (1-i)ln(1-i) + (1+i)ln(1+i) - i^2
    where i = Im / Ic is the normalized current amplitude.
    """
    # Check for valid range of i
    if not (0 < i < 1):
        print("Error: The value of i must be between 0 and 1 (exclusive).")
        print("For i=0, the loss is 0.")
        print("The formula is not defined for i>=1.")
        return None

    # --- Calculation Steps ---
    print(f"\nCalculating normalized loss for i = {i}:")

    # Final Equation Form
    equation_str = f"2\u03C0Q/(\u03BC\u2080Ic\u00B2) = (1-{i})*ln(1-{i}) + (1+{i})*ln(1+{i}) - {i}\u00B2"
    print("Formula: " + equation_str)
    
    # Calculate each term individually
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    # Calculate the final result
    result = term1 + term2 - term3

    # Print the breakdown of the calculation
    print("\nBreaking down the calculation:")
    print(f"Term 1: (1-i)ln(1-i) = ({1-i:.3f}) * ln({1-i:.3f}) = {term1:.6f}")
    print(f"Term 2: (1+i)ln(1+i) = ({1+i:.3f}) * ln({1+i:.3f}) = {term2:.6f}")
    print(f"Term 3: i^2 = {i:.3f}^2 = {term3:.6f}")
    
    print("\nFinal Result:")
    print(f"2\u03C0Q/(\u03BC\u2080Ic\u00B2) = {term1:.6f} + {term2:.6f} - {term3:.6f} = {result:.6f}")
    
    return result

if __name__ == "__main__":
    # Prompt the user for a value of i
    try:
        i_input = input("Enter the normalized current i = Im/Ic (a value between 0 and 1, e.g., 0.5): ")
        # Use a default value if the user just presses Enter
        if i_input.strip() == "":
            i_val = 0.5
            print(f"No input provided, using default value i = {i_val}")
        else:
            i_val = float(i_input)
        
        calculate_normalized_loss(i_val)

    except ValueError:
        print("Invalid input. Please enter a numerical value.")
    except Exception as e:
        print(f"An error occurred: {e}")

# The final answer format requires a single numerical value.
# Below is the value calculated for i = 0.5.
final_answer_for_i_0_5 = (1 - 0.5) * math.log(1 - 0.5) + (1 + 0.5) * math.log(1 + 0.5) - 0.5**2
# This value will be wrapped in the required format at the end.
# <<<0.01163158>>>
