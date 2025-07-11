import math

def calculate_normalized_ac_loss(i):
    """
    Calculates the normalized AC loss for a superconducting elliptic bar
    with a transport AC current for i = Im/Ic < 1.

    The function calculates the value of 2 * pi * Q / (mu_0 * Ic^2), where Q is
    the loss per cycle per unit length.

    The formula used is from W. T. Norris:
    Loss = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]
    
    Args:
        i (float): The ratio of the current amplitude to the critical current, Im/Ic.
                   Must be between 0 (exclusive) and 1 (exclusive).
    """
    # Validate the input parameter i
    if not (0 < i < 1):
        print("Error: The value of i must be between 0 and 1 (0 < i < 1).")
        return

    # The equation for the normalized loss: 2 * [ (1-i)ln(1-i) + (1+i)ln(1+i) - i^2 ]
    # We will break it down and print each component.

    print(f"Calculating the normalized loss for i = {i}:")
    print("Equation: 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]")
    print("-" * 50)

    # First term: (1-i)*ln(1-i)
    val_1_minus_i = 1 - i
    log_val_1_minus_i = math.log(val_1_minus_i)
    term1 = val_1_minus_i * log_val_1_minus_i
    print(f"Term (1-i)*ln(1-i):")
    print(f"  1 - i = {val_1_minus_i:.4f}")
    print(f"  ln(1-i) = {log_val_1_minus_i:.4f}")
    print(f"  Value = ({val_1_minus_i:.4f}) * ({log_val_1_minus_i:.4f}) = {term1:.4f}")
    print("-" * 50)
    
    # Second term: (1+i)*ln(1+i)
    val_1_plus_i = 1 + i
    log_val_1_plus_i = math.log(val_1_plus_i)
    term2 = val_1_plus_i * log_val_1_plus_i
    print(f"Term (1+i)*ln(1+i):")
    print(f"  1 + i = {val_1_plus_i:.4f}")
    print(f"  ln(1+i) = {log_val_1_plus_i:.4f}")
    print(f"  Value = ({val_1_plus_i:.4f}) * ({log_val_1_plus_i:.4f}) = {term2:.4f}")
    print("-" * 50)

    # Third term: -i^2
    term3 = i**2
    print(f"Term -i^2:")
    print(f"  i^2 = {i:.4f}^2 = {term3:.4f}")
    print(f"  Value = -{term3:.4f}")
    print("-" * 50)

    # Sum of terms inside the bracket
    bracket_sum = term1 + term2 - term3
    print("Combining terms inside the bracket:")
    print(f"  ({term1:.4f}) + ({term2:.4f}) - ({term3:.4f}) = {bracket_sum:.4f}")
    print("-" * 50)

    # Final result
    final_result = 2 * bracket_sum
    print("Final Normalized Loss Result [2*pi*Q / (mu_0*Ic^2)]:")
    print(f"  2 * ({bracket_sum:.4f}) = {final_result:.6f}")
    
    return final_result

if __name__ == '__main__':
    # Set an example value for the normalized current i
    # You can change this value to any number between 0 and 1.
    normalized_current = 0.5
    
    # Calculate and print the normalized loss
    result = calculate_normalized_ac_loss(normalized_current)
    
    # Final answer in the specified format
    if result is not None:
        print(f"\n<<<answer content>>>\n{result}")
