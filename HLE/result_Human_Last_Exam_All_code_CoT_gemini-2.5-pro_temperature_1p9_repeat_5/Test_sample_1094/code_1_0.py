import math

def calculate_normalized_ac_loss(i):
    """
    Calculates the normalized AC loss for an elliptical superconductor.

    The function calculates the loss per cycle Q for a transport AC current
    of amplitude Im, in the normalized form 2*pi*Q / (mu_0 * Ic^2), where
    i = Im / Ic < 1 is the normalized current.

    The calculation is based on the Norris formula for an elliptical bar.
    """

    # --- Input Validation ---
    if not (0 <= i < 1):
        print("Error: The normalized current 'i' must be in the range [0, 1).")
        print("The formula is not valid for i >= 1.")
        return

    # --- Introduction ---
    print(f"Calculating normalized AC loss for a normalized current i = {i}\n")
    print("The normalized loss, Q_norm = 2*pi*Q / (mu_0 * Ic^2), is given by the formula:")
    print("Q_norm = 2 * [ (1 - i)ln(1 - i) + (1 + i)ln(1 + i) - i^2 ]\n")

    # --- Calculation of Each Term ---
    print("Step-by-step calculation:")

    # Term 1: (1 - i) * ln(1 - i)
    # Note: math.log is the natural logarithm (ln)
    term1_factor1 = 1 - i
    term1_factor2 = math.log(1 - i)
    term1 = term1_factor1 * term1_factor2
    print(f"Term 1: (1 - {i}) * ln(1 - {i}) = {term1_factor1:.4f} * {term1_factor2:.4f} = {term1:.6f}")

    # Term 2: (1 + i) * ln(1 + i)
    term2_factor1 = 1 + i
    term2_factor2 = math.log(1 + i)
    term2 = term2_factor1 * term2_factor2
    print(f"Term 2: (1 + {i}) * ln(1 + {i}) = {term2_factor1:.4f} * {term2_factor2:.4f} = {term2:.6f}")

    # Term 3: -i^2
    term3 = i**2
    print(f"Term 3: i^2 = {i}^2 = {term3:.6f}\n")

    # --- Final Calculation ---
    # Sum of terms inside the bracket
    bracket_sum = term1 + term2 - term3
    print("Sum of terms inside the main bracket:")
    print(f"({term1:.6f}) + ({term2:.6f}) - ({term3:.6f}) = {bracket_sum:.6f}\n")

    # Final result
    normalized_loss = 2 * bracket_sum
    print("Final Result:")
    print(f"Q_norm = 2 * {bracket_sum:.6f} = {normalized_loss:.6f}")
    
    # Returning the final numerical answer for the program
    return normalized_loss


if __name__ == '__main__':
    # Set the value for the normalized current i = Im/Ic.
    # You can change this value to any number between 0 and 1.
    normalized_current_i = 0.7

    # Calculate and print the result.
    final_answer = calculate_normalized_ac_loss(normalized_current_i)

    # The final answer in the required format
    # The string representation is formatted to avoid floating point inaccuracies
    # that might arise in different environments.
    if final_answer is not None:
        print(f"\n<<<{final_answer:.6f}>>>")
