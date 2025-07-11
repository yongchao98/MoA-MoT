import decimal

def solve_tan_googol():
    """
    Calculates the first 3 digits after the comma of tan(10^100).
    """
    # 1. Set a high precision for the decimal calculations.
    # 10^100 has 101 digits. The number of integer multiples of pi in it
    # is about 100 digits long. We need precision greater than that.
    # 115 digits is a safe level of precision for this task.
    decimal.getcontext().prec = 115

    # 2. Define the number N = 10^100 and a high-precision value of pi.
    N = decimal.Decimal(10) ** 100
    
    # A high-precision string representation of pi.
    pi_str = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
    pi = decimal.Decimal(pi_str)

    # 3. Calculate the remainder of N divided by pi.
    # This leverages the periodicity of tan(x), where tan(x) = tan(x mod pi).
    remainder = N % pi

    # 4. Calculate the tangent of the remainder.
    result = remainder.tan()

    # 5. Output the numbers in the final equation as requested.
    print(f"tan(10^100) is equivalent to tan(10^100 mod pi)")
    print(f"pi ≈ {pi}")
    print(f"10^100 mod pi = {remainder}")
    print(f"tan(10^100) ≈ {result}")

    # 6. Extract the first 3 digits after the decimal point.
    result_str = f"{result:.10f}" # Format to ensure enough decimal places
    decimal_point_index = result_str.find('.')
    first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

    print("\n---")
    print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")


solve_tan_googol()