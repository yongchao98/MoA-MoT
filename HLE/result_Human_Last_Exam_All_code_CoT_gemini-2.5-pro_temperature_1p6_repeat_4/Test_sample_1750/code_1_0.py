import math

def calculate_integral_fraction():
    """
    This function calculates the fractional part of the integral result.
    The integral I simplifies to (pi / 2**49) * C(50, 25).
    This script calculates the numerator and denominator of the rational
    multiplier of pi.
    """
    
    # Calculate C(50, 25)
    # C(n, k) = n! / (k! * (n-k)!)
    n = 50
    k = 25
    
    # math.comb is available from Python 3.8
    try:
        numerator = math.comb(n, k)
    except AttributeError:
        # Manual implementation for older Python versions
        if k < 0 or k > n:
            numerator = 0
        if k == 0 or k == n:
            numerator = 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        numerator = res

    # Denominator is 2^49
    denominator = 2**49
    
    # We can simplify the fraction by finding the greatest common divisor.
    # Numerator C(50, 25) is an integer. Let's find how many factors of 2 it has.
    # Legendre's formula: v_2(n!) = n - s_2(n) where s_2(n) is sum of binary digits.
    # v_2(50!) = 50 - (1+1+0+0+1+0) = 45
    # v_2(25!) = 25 - (1+1+0+0+1) = 22
    # v_2(C(50,25)) = v_2(50!) - 2*v_2(25!) = 45 - 2*22 = 1.
    # So C(50, 25) has one factor of 2.
    
    simplified_numerator = numerator // 2
    simplified_denominator = denominator // 2
    
    print("The integral evaluates to (N/D) * pi")
    print(f"The numerator N is: {numerator}")
    print(f"The denominator D is: 2**49 = {denominator}")
    print(f"After simplification, the fraction is ({simplified_numerator}/{simplified_denominator}) * pi")
    print("Final answer in the requested format will be the fraction part.")
    
    # The final equation is I = (numerator / denominator) * pi
    # The problem asks to output each number in the final equation.
    # Let's consider the final answer to be the fraction multiplier of pi.
    final_answer_fraction = f"{simplified_numerator}/{simplified_denominator}"
    # The problem asks for the answer directly at the end.
    # Let's provide the simplified fraction string as the answer.

calculate_integral_fraction()