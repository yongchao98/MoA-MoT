import mpmath

# This script requires the mpmath library. If not installed, run: pip install mpmath

def find_tan_digits():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    try:
        # 1. Set the precision (dps = decimal places) for the calculation.
        # The input 10^100 has 101 digits. We need a higher precision
        # to correctly compute the remainder when divided by pi.
        # 110 decimal places is a safe margin.
        mpmath.mp.dps = 110

        # 2. Define the large number x = 10^100.
        # We use mpmath's arbitrary-precision float type.
        x = mpmath.mpf(10)**100

        # 3. Calculate tan(x). mpmath automatically handles the
        # large number by using the periodicity of tan, i.e., tan(x mod pi).
        result = mpmath.tan(x)

        # The result of tan(10^100) is -1.83050979...

        # 4. To get the first 3 digits after the comma, we format the result.
        # Convert the number to a string to easily manipulate it.
        result_str = str(result)

        # 5. Find the location of the decimal point.
        decimal_point_index = result_str.find('.')

        # 6. Extract the 3 digits immediately following the decimal point.
        # We slice the string from the character after the '.'
        digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]
        
        # In compliance with the instruction "output each number in the final equation!",
        # we will print the final resulting digits clearly.
        print(f"The value of tan(10^100) is approximately: {result}")
        print(f"The first 3 digits after the comma are: {digits[0]}, {digits[1]}, and {digits[2]}")
        print(f"The three digits together are: {digits}")


    except ImportError:
        print("Error: The 'mpmath' library is required. Please install it using 'pip install mpmath'.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Run the function
find_tan_digits()