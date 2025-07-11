import decimal

def find_tan_googol_digits():
    """
    Calculates the first 3 digits after the decimal point of tan(10^100).
    """
    # Set a high precision for the decimal calculations.
    # To find `10^100 mod pi` accurately, we need pi to be precise
    # to more than 100 decimal places. 120 should be sufficient.
    prec = 120
    decimal.getcontext().prec = prec

    # A high-precision string of Pi with 119 decimal places (120 significant digits).
    pi_str = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"

    # --- High-precision math functions for Decimal objects ---

    def sin_decimal(x, p):
        """Computes sin(x) for a Decimal 'x' to 'p' places using Taylor series."""
        decimal.getcontext().prec = p + 4  # Use a few guard digits
        x_squared = x * x
        result = decimal.Decimal(0)
        term = x
        i = 1
        # Taylor series for sin(x) = x - x^3/3! + x^5/5! - ...
        while abs(term) > decimal.Decimal('1e-' + str(p + 2)):
            result += term
            i += 2
            term *= -x_squared / ((i - 1) * i)
        decimal.getcontext().prec = p
        return +result

    def cos_decimal(x, p):
        """Computes cos(x) for a Decimal 'x' to 'p' places using Taylor series."""
        decimal.getcontext().prec = p + 4  # Use a few guard digits
        x_squared = x * x
        result = decimal.Decimal(0)
        term = decimal.Decimal(1)
        i = 0
        # Taylor series for cos(x) = 1 - x^2/2! + x^4/4! - ...
        while abs(term) > decimal.Decimal('1e-' + str(p + 2)):
            result += term
            i += 2
            term *= -x_squared / ((i - 1) * i)
        decimal.getcontext().prec = p
        return +result

    def tan_decimal(x, p):
        """Computes tan(x) for a Decimal 'x' to 'p' places."""
        s = sin_decimal(x, p + 4)
        c = cos_decimal(x, p + 4)
        if c == 0:
            raise ValueError("Tangent is undefined.")
        decimal.getcontext().prec = p
        return s / c

    # --- Main calculation ---

    # 1. Initialize Decimal objects for Pi and the googol (10^100).
    pi = decimal.Decimal(pi_str)
    # The final equation includes 10 and 100.
    power_base = decimal.Decimal(10)
    exponent = decimal.Decimal(100)
    googol = power_base ** exponent

    # 2. Calculate r = 10^100 mod pi. The result is in [0, pi).
    r = googol % pi

    # 3. Reduce r to the interval [-pi/2, pi/2] for faster convergence of the
    #    Taylor series, using the identity tan(x) = tan(x - pi).
    pi_half = pi / 2
    if r > pi_half:
        r = r - pi

    # 4. Calculate tan(r). We need 3 digits, so 10-digit precision is more than enough.
    tan_prec = 10
    final_tan_value = tan_decimal(r, tan_prec)

    # 5. Format the result and extract the required digits.
    # We want the first 3 digits *after* the decimal point.
    tan_str = format(final_tan_value, f'.{tan_prec}f')
    
    decimal_point_index = tan_str.find('.')
    if decimal_point_index != -1:
        first_three_digits = tan_str[decimal_point_index + 1 : decimal_point_index + 4]
        # Printing the final equation with its components and the answer
        print(f"The value of tan({int(power_base)}^{int(exponent)}) is approximately: {final_tan_value}")
        print(f"The first 3 digits after the decimal point are: {first_three_digits}")
    else:
        print("Could not find decimal point in the result string.")


# Execute the function to find and print the result.
find_tan_googol_digits()