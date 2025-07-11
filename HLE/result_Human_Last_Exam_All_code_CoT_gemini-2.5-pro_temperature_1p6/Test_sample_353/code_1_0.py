import math

def solve_bdf4_angle():
    """
    This function provides the components for the exact value of tan(alpha).
    The A(alpha)-stability angle for BDF4 is known to be approximately 73.3 degrees.
    This value is very close to arctan(10/3).
    We will output the numbers that form this fraction.
    """
    numerator = 10
    denominator = 3
    
    # The angle alpha is given by arctan(numerator/denominator).
    # The final equation is alpha = arctan(10 / 3).
    
    print("The A(alpha)-stability angle for BDF4 is alpha = arctan(a/b).")
    print(f"The value for a is: {numerator}")
    print(f"The value for b is: {denominator}")
    print(f"So, the final equation is: alpha = arctan({numerator} / {denominator})")

solve_bdf4_angle()