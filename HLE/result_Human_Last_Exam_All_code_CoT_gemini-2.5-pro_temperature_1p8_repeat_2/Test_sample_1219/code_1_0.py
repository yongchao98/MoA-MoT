import decimal

def solve_tan_gogol():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    try:
        # 1. Set a high precision for the calculation.
        # 10^100 has 101 digits, so we need more than that.
        # Let's set the precision to 110 digits for accuracy.
        decimal.getcontext().prec = 110

        # 2. Define a high-precision string for Pi. This is necessary because
        # the result of 10^100 mod pi depends on the precision of pi.
        pi_str = "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647"
        pi = decimal.Decimal(pi_str)
        
        # 3. The number we want to find the tangent of, as a Decimal object.
        x = decimal.Decimal(10) ** 100

        # 4. Use the periodicity of tan: tan(x) = tan(x mod pi).
        # We calculate the remainder of 10^100 divided by pi.
        remainder = x % pi

        # 5. Calculate the tangent of the remainder.
        tan_value = remainder.tan()
        
        # As per the instruction, we show the numbers in the final equation.
        # Here, the conceptual equation is tan(10^100) = value
        print(f"The full result for tan(10^100) is: {tan_value}")
        
        # 6. Extract the first 3 digits after the comma.
        # Convert the Decimal result to a string.
        tan_string = str(tan_value)
        
        # Find the position of the decimal point.
        decimal_point_index = tan_string.find('.')
        
        # Slice the string to get the 3 digits right after the decimal point.
        if decimal_point_index != -1:
            start_index = decimal_point_index + 1
            three_digits = tan_string[start_index : start_index + 3]
            print(f"The first 3 digits after the comma of tan(10^100) are: {three_digits}")
        else:
            # Handle cases where there might not be a decimal part, which is unlikely here.
            print("The result is an integer, so there are no digits after the comma.")

    except Exception as e:
        print(f"An error occurred: {e}")

solve_tan_gogol()