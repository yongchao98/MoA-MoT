import sys
from decimal import Decimal, getcontext

def solve_tan_large_number():
    """
    Calculates the first 3 digits after the comma of tan(10^100)
    using high-precision arithmetic.
    """
    # The argument 10^100 has 101 digits. To find the remainder when dividing
    # by pi, we need precision greater than 101. Let's set it to 120 for safety.
    try:
        getcontext().prec = 120
    except ValueError:
        print("Fatal error: Could not set the required decimal precision (120).")
        print("Your Python environment may have a lower limit.")
        sys.exit(1)

    # We need to compute tan(10^100).
    # Since tan(x) has a period of pi, this is equal to tan(10^100 mod pi).
    # We will use the decimal module to perform this calculation with high precision.
    
    # The numbers in the final equation tan(10**100) = tan(rem) = result
    
    # 1. The original number's exponent
    power = 100
    
    # 2. The high-precision representation of 10^100
    x = Decimal(10) ** power
    
    # 3. A high-precision value for pi
    pi = 4 * Decimal(1).atan()
    
    # 4. The argument for tan(), reduced modulo pi
    reduced_argument = x % pi
    
    # 5. The final result of the tangent function
    result = reduced_argument.tan()

    print(f"The calculation follows the identity: tan(10**{power}) = tan(10**{power} % pi)")
    print(f"\nThe reduced argument (10**{power} % pi) is: {reduced_argument}")
    print(f"The final result tan(reduced_argument) is: {result}")
    
    # Convert the result to a string to find the digits after the decimal point.
    # We format it to a sufficient number of decimal places to avoid scientific notation.
    result_str = f"{result:.50f}"
    
    # Find the position of the decimal point
    if '.' in result_str:
        decimal_point_index = result_str.find('.')
        # Extract the first three digits after the decimal point
        first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]
    else:
        # This case is unlikely for a tan result but is handled for completeness.
        first_three_digits = "000"

    print(f"\nThus, the first 3 digits after the comma of tan(10^100) are: {first_three_digits}")

solve_tan_large_number()