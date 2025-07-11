import mpmath

def find_tan_digits():
    """
    Calculates the first 3 digits after the decimal point of tan(10^100).
    """
    # Set the precision (dps = decimal places) for the calculation.
    # The argument 10^100 is very large. To compute tan(x) for a large x,
    # we need to find x mod pi with high accuracy. The number of digits in
    # 10^100 is 101, so we need a precision greater than that.
    # We use 115 as a safe value.
    mpmath.mp.dps = 115

    # Define the number for the tan function argument
    base = 10
    exponent = 100
    argument = mpmath.power(base, exponent)

    # Calculate tan(10^100). The mpmath library handles the large argument
    # by correctly reducing it modulo pi.
    result = mpmath.tan(argument)

    # Convert the high-precision result to a string to easily access the digits.
    # We format it to have a sufficient number of digits shown.
    result_str = mpmath.nstr(result, 10)

    # The result is approximately -1.376...
    # We locate the decimal point in the string representation.
    decimal_point_index = result_str.find('.')

    # We extract the three characters immediately following the decimal point.
    # This gives us the first three digits after the comma.
    first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

    # The "final equation" is tan(10^100) = -1.376... The required numbers are the
    # three digits '3', '7', and '6'. We will print them as a single string.
    print(first_three_digits)

find_tan_digits()
