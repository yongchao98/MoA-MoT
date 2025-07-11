from fractions import Fraction
from math import comb

def double_factorial(n):
    """Computes the double factorial n!!."""
    if not isinstance(n, int):
        raise TypeError("Input must be an integer.")
    if n < -1:
        raise ValueError("Double factorial is not defined for integers less than -1.")
    if n == -1 or n == 0:
        return 1
    
    res = 1
    i = n
    while i >= 1:
        res *= i
        i -= 2
    return res

def solve_integral():
    """
    This function calculates the rational coefficient of pi for the integral.
    The integral is given by the formula:
    I = pi * sum_{j=0 to 25} C(50, 2j) * 2^(50-2j) * (49!! * (49-2j)!!) / (100-2j)!!
    """
    n = 25
    total_sum = Fraction(0)
    
    print("The integral is equal to pi times the sum of the following terms:")
    
    # Pre-calculate the constant part of the numerator
    df_num1 = double_factorial(2 * n - 1) # 49!!
    
    for j in range(n + 1):
        # Binomial coefficient C(50, 2j)
        comb_term = comb(2 * n, 2 * j)
        
        # Power of 2 term: 2^(50-2j)
        pow_term = 2**(2 * n - 2 * j)
        
        # Double factorial terms
        df_num2 = double_factorial(2 * n - 2 * j - 1) # (49-2j)!!
        df_den = double_factorial(4 * n - 2 * j)     # (100-2j)!!
        
        # Calculate the numerator and denominator of the term
        numerator = comb_term * pow_term * df_num1 * df_num2
        denominator = df_den
        
        term = Fraction(numerator, denominator)
        
        # Print each term in the final equation
        if j > 0:
            print(" + ", end="")
        print(f"{term}", end="")

        total_sum += term
        
    print(f"\n\nThe sum evaluates to: {total_sum}")
    print(f"\nTherefore, the integral is: {total_sum.numerator}/{total_sum.denominator} * pi")
    
    # Final answer as per the requested format
    print("\n<<<", end="")
    print(f"{total_sum.numerator}/{total_sum.denominator}", end="")
    print(">>>")


if __name__ == '__main__':
    solve_integral()
