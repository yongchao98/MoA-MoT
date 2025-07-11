import math

def calculate_pdf_value(z):
    """
    Calculates the value of the probability density function f_Z(z) at a given point z.

    The PDF is derived from a problem involving order statistics and is given by
    f_Z(z) = 30 * z^2 * (1-z)^2.
    """
    # The coefficient in the PDF formula
    coefficient = 30

    # Calculate z^2
    z_squared = z**2

    # Calculate (1-z)
    one_minus_z = 1 - z

    # Calculate (1-z)^2
    one_minus_z_squared = one_minus_z**2

    # Calculate the final PDF value
    f_z = coefficient * z_squared * one_minus_z_squared

    # Print the step-by-step calculation
    print(f"The task is to calculate f_Z({z}) for the given probability density function.")
    print(f"The formula for the PDF is: f_Z(z) = {coefficient} * z^2 * (1-z)^2")
    print(f"The final numerical equation is: {coefficient} * {z}^2 * (1-{z})^2 = Result")
    print("\nCalculation steps:")
    print(f"1. First term: {coefficient}")
    print(f"2. Second term: {z}^2 = {z_squared}")
    print(f"3. Third term: (1-{z})^2 = {one_minus_z}^2 = {round(one_minus_z_squared, 4)}") # Round for cleaner output
    print(f"4. The final equation is: {coefficient} * {z_squared} * {round(one_minus_z_squared, 4)} = {f_z}")

    # The user prompt asks for print, not return, but for clarity we show the final answer
    print("\nFinal Answer:")
    print(f"f_Z({z}) = {f_z}")

# The specific value we need to calculate
z_value = 0.2
calculate_pdf_value(z_value)