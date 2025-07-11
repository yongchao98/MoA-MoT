import decimal

def sin(x):
    """
    Computes sin(x) for a Decimal object x using its Taylor series.
    sin(x) = x - x^3/3! + x^5/5! - ...
    """
    # Temporarily increase precision for intermediate calculations
    decimal.getcontext().prec += 2
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    # Series evaluation continues until the sum stops changing
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    # Restore original precision
    decimal.getcontext().prec -= 2
    return +s

def cos(x):
    """
    Computes cos(x) for a Decimal object x using its Taylor series.
    cos(x) = 1 - x^2/2! + x^4/4! - ...
    """
    decimal.getcontext().prec += 2
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    decimal.getcontext().prec -= 2
    return +s

def tan(x):
    """
    Computes tan(x) for a Decimal object x.
    tan(x) = sin(x) / cos(x)
    """
    return sin(x) / cos(x)

def find_tan_googol_digits():
    """
    Calculates tan(10^100) and finds the first 3 digits after the decimal point.
    """
    try:
        # 1. Set precision. 10^100 has 101 digits, so we need more to find the remainder
        # accurately. 120 digits provides a safe margin.
        decimal.getcontext().prec = 120

        # 2. Define a high-precision value of Pi as a string.
        # This string must have at least as many digits as the precision.
        pi_str = "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647"
        pi = decimal.Decimal(pi_str)

        # 3. Represent 10^100 as a Decimal object.
        x = decimal.Decimal(10) ** 100

        # 4. Find 10^100 mod pi.
        # This is the remainder r = x - n*pi, where n = floor(x/pi).
        # This remainder will be the argument to our tan function.
        q = (x / pi).to_integral_value(rounding=decimal.ROUND_FLOOR)
        argument = x - q * pi

        # 5. Calculate tan(argument).
        result = tan(argument)

        # 6. Output the full equation and result as requested.
        print(f"The calculation is for tan(10^100).")
        print(f"The argument 10^100 is reduced to {argument} using periodicity.")
        print(f"tan(10^100) = {result}")
        
        # 7. Extract the first 3 digits after the comma.
        s_result = str(result)
        if '.' in s_result:
            fractional_part = s_result.split('.')[1]
            first_three_digits = (fractional_part + "000")[:3]
        else:
            first_three_digits = "000"

        print(f"\nThe first 3 digits after the comma of tan(10^100) are: {first_three_digits}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_tan_googol_digits()
