import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for a superconducting elliptic bar.

    The function calculates the loss per cycle Q, normalized as 2*pi*Q/(mu_0*Ic^2),
    for a transport AC current with amplitude Im.

    Args:
        i (float): The ratio of the current amplitude to the critical current, i = Im/Ic.
                   Must be in the range [0, 1).

    Returns:
        float: The value of the normalized loss, or None if i is out of range.
    """
    if not 0 <= i < 1:
        print("Error: The value of i must be in the range [0, 1).")
        return None

    # Handle the edge case i=0 to avoid potential floating point issues, though the formula works.
    if i == 0:
        return 0.0

    # The normalized loss is given by the function:
    # f(i) = 2 * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i^2 ]

    # Breaking down the calculation to show each number
    print("Calculating normalized loss for i = {}".format(i))
    print("-" * 35)

    # First term: 2 * (1 - i) * ln(1 - i)
    term1_factor = 1 - i
    term1_log = math.log(term1_factor)
    term1 = term1_factor * term1_log
    print("Term (1 - i):              {:.6f}".format(term1_factor))
    print("Term ln(1 - i):            {:.6f}".format(term1_log))
    print("Product (1-i)*ln(1-i):     {:.6f}".format(term1))
    print("-" * 35)

    # Second term: 2 * (1 + i) * ln(1 + i)
    term2_factor = 1 + i
    term2_log = math.log(term2_factor)
    term2 = term2_factor * term2_log
    print("Term (1 + i):              {:.6f}".format(term2_factor))
    print("Term ln(1 + i):            {:.6f}".format(term2_log))
    print("Product (1+i)*ln(1+i):     {:.6f}".format(term2))
    print("-" * 35)

    # Third term: -2 * i^2
    term3_square = i**2
    print("Term i^2:                  {:.6f}".format(term3_square))
    print("-" * 35)

    # Combine terms inside the bracket
    bracket_value = term1 + term2 - term3_square
    print("Value inside brackets [...]: {:.6f}".format(bracket_value))

    # Final result
    result = 2 * bracket_value
    return result

if __name__ == "__main__":
    # Example usage: Calculate the normalized loss for i = 0.6
    # You can change this value to any float between 0 and 1.
    i_value = 0.6

    print("The normalized AC loss 2\u03C0Q/(\u03BC\u2080I\u209C\u00B2) for an elliptic bar is given by the function:")
    print("f(i) = 2 * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i\u00B2 ] for i < 1\n")

    normalized_loss = calculate_normalized_loss(i_value)

    if normalized_loss is not None:
        print("\n=======================================================")
        print("Final Result: 2\u03C0Q/(\u03BC\u2080I\u209C\u00B2) = {:.6f}".format(normalized_loss))
        print("=======================================================")
