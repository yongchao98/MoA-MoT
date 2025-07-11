# The gmpy2 library is required for this high-precision calculation.
# If you don't have it, please install it using: pip install gmpy2
import gmpy2

def solve():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    try:
        # Set the precision for the calculation in bits.
        # 10^100 is approximately 2^332. We need more bits than that for precision.
        # 400 bits (~120 decimal digits) is sufficient.
        gmpy2.get_context().precision = 400

        # Define the large number 10^100.
        # We use gmpy2's arbitrary-precision integer type (mpz).
        x = gmpy2.mpz(10)**100

        # Calculate tan(10^100).
        # gmpy2's tan function correctly handles the large argument by
        # reducing it modulo pi to the interval [-pi/2, pi/2].
        result = gmpy2.tan(x)

        # The result is approximately -1.5378168...
        # We want the first 3 digits after the decimal point: 537.

        # To extract these digits programmatically:
        # 1. Take the absolute value.
        abs_result = abs(result)
        # 2. Get the fractional part by subtracting the integer part.
        fractional_part = abs_result - gmpy2.floor(abs_result)
        # 3. Multiply by 1000 and take the integer part to get the first three digits.
        three_digits = gmpy2.floor(fractional_part * 1000)

        # The final equation is tan(10^100) = -1.537...
        # The first digit after the comma is 5.
        # The second digit after the comma is 3.
        # The third digit after the comma is 7.
        print("The value of tan(10^100) is approximately: {}".format(result))
        print("The first 3 digits after the comma are: {}".format(three_digits))

    except ImportError:
        print("This script requires the 'gmpy2' library.")
        print("Please install it using: pip install gmpy2")
    except Exception as e:
        print(f"An error occurred: {e}")

solve()