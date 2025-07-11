import decimal

def solve_tan_large_number():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # Step 1: Set the precision for the decimal calculations.
    # 10^100 is a number with 101 digits. To accurately compute its remainder
    # when divided by pi, our working precision must be greater than 101.
    # We'll set it to 110 to have enough guard digits for the calculation.
    decimal.getcontext().prec = 110

    # Step 2: Use a high-precision value of pi.
    # This string contains pi to more than 110 digits.
    pi = decimal.Decimal("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679")

    # Step 3: Represent 10^100 as a Decimal object.
    x = decimal.Decimal(10)**100

    # Step 4: Use the periodicity of tan(x) = tan(x mod pi)
    # We calculate y = 10^100 mod pi.
    y = x % pi

    # Step 5: Calculate the tangent of the remainder.
    result = y.tan()

    # Step 6: Extract the first 3 digits after the comma.
    # The result will be a string like '-1.486...'. We split it at the
    # decimal point and take the first three characters of the fractional part.
    fractional_part = str(result).split('.')[1]
    first_three_digits = fractional_part[:3]

    # As requested, output the numbers in the final equation.
    print(f"The equation we want to solve is: tan(10^100)")
    print(f"Using the periodicity of the tangent function, this is equal to: tan(10^100 mod pi)")
    print(f"Value of pi used (first 50 digits): {str(pi)[:52]}...")
    print(f"Value of (10^100 mod pi): {y}")
    print(f"Result of tan(10^100): {result}")
    print(f"The first 3 digits after the comma are: {first_three_digits}")

solve_tan_large_number()