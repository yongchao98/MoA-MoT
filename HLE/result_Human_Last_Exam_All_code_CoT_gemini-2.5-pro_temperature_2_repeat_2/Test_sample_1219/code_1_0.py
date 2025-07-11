import mpmath

# Before running, you might need to install the mpmath library:
# pip install mpmath

def solve_tan_googol():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # Set the working precision for the calculation.
    # To compute tan(10^100), we need to find 10^100 mod pi.
    # This requires pi to be known to a precision greater than the number of digits in 10^100.
    # 10^100 has 101 digits. Let's use a precision of 120 decimal places for safety.
    mpmath.mp.dps = 120

    # The input number, a googol
    x = 10**100

    # Calculate the tangent. mpmath automatically handles the large argument by
    # reducing it modulo pi.
    result = mpmath.tan(x)

    # Convert the high-precision result to a string with plenty of digits.
    result_str = mpmath.nstr(result, 10)

    # Split the string at the decimal point to isolate the fractional part.
    if '.' in result_str:
        decimal_part = result_str.split('.')[1]
        # The first three digits after the comma are the first three characters of the decimal part.
        first_three_digits = decimal_part[:3]
    else:
        # This case handles an integer result, where the decimal digits would be '000'.
        first_three_digits = '000'

    # Print the final equation with the computed value.
    print(f"tan(10^100) = {result_str}")
    
    # Print the requested first three digits.
    print(f"The first 3 digits after the comma are: {first_three_digits}")

if __name__ == "__main__":
    solve_tan_googol()
