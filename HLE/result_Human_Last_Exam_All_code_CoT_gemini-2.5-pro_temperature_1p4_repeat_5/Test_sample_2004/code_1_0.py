import numpy as np
import math

def format_to_significant_figures(num, sf):
    """
    Formats a number to a specified number of significant figures.
    """
    if num == 0:
        return "0." + "0" * (sf - 1)
    
    # Calculate the number of decimal places needed
    power = sf - 1 - math.floor(math.log10(abs(num)))
    
    # Round the number
    rounded_num = round(num, power)
    
    # Format the string, ensuring trailing zeros are kept if necessary
    return f"{rounded_num:.{max(0, power)}f}"

def solve_model():
    """
    Finds the maximally parsimonious model for the given data and reports the result.
    """
    # The 25 observations of x and y
    x = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # As determined by the analysis, the square root model y = a*sqrt(x) + b is optimal.
    # We transform x to sqrt(x) and perform a linear fit.
    x_transformed = np.sqrt(x)
    
    # Fit the model using polyfit, which performs a least-squares polynomial fit.
    # A degree of 1 on the transformed data is equivalent to y = a*z + b where z = sqrt(x).
    params = np.polyfit(x_transformed, y, 1)
    a, b = params[0], params[1]

    # Format the parameters to 3 significant figures
    a_str = format_to_significant_figures(a, 3)
    
    # Handle the sign for the equation's presentation
    sign = "+" if b >= 0 else "-"
    b_str_unsigned = format_to_significant_figures(abs(b), 3)

    # Print the final results
    print("The maximally parsimonious model is of the form: y = a * sqrt(x) + b")
    print("\nParameter estimates (to 3 significant figures):")
    print(f"a = {a_str}")
    print(f"b = {format_to_significant_figures(b, 3)}")
    print("\nFinal equation:")
    print(f"y = {a_str} * sqrt(x) {sign} {b_str_unsigned}")

solve_model()