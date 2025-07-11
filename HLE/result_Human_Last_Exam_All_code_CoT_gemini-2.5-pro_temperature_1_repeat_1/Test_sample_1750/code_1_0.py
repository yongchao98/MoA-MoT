import math

def calculate_integral_fraction():
    """
    This function calculates the rational part of the integral's result based on a known formula.
    
    The integral simplifies to:
    I = integral from 0 to pi of (sin(4x) - sin(2x))^50 dx

    A known formula for this type of integral (for even power n=50) is:
    I = (pi / 2**(n-1)) * C(n, n/2)
    where C(n, k) is the binomial coefficient "n choose k".

    For n = 50, the result is (pi / 2**49) * C(50, 25).
    The user wants the answer as a fraction, so we will compute the rational part C(50, 25) / 2**49.
    """
    n = 50
    k = n // 2

    # Calculate the binomial coefficient C(n, k) = C(50, 25)
    # math.comb is available in Python 3.8+
    try:
        numerator = math.comb(n, k)
    except AttributeError:
        # Manual implementation for older Python versions
        if k < 0 or k > n:
            numerator = 0
        elif k == 0 or k == n:
            numerator = 1
        elif k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        numerator = res
        
    # Calculate the denominator 2**(n-1) = 2**49
    denominator = 2**(n - 1)

    # The problem asks to show the numbers in the final equation.
    print("The integral simplifies to the form: (pi / 2^(n-1)) * C(n, n/2)")
    print(f"For n=50, this is: (pi / 2^49) * C(50, 25)")
    print("\nThe equation for the integral's value is:")
    print(f"I = ({numerator} / {denominator}) * pi")
    
    # The problem asks for the answer as a strict fraction. This implies returning the rational part.
    # We find the greatest common divisor (GCD) to simplify the fraction.
    common_divisor = math.gcd(numerator, denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    print("\nThe rational part of the answer, as a simplified fraction, is:")
    print(f"{simplified_numerator}/{simplified_denominator}")

calculate_integral_fraction()