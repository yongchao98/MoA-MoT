import math
from fractions import Fraction

def double_factorial(n):
    """Computes the double factorial of n."""
    if n < 0:
        return 1 # By convention for this problem's formula, (-1)!! = 1
    if n == 0:
        return 1
    res = 1
    for i in range(n, 0, -2):
        res *= i
    return res

def combinations(n, k):
    """Computes n choose k."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def main():
    """
    Calculates the sum which is the rational coefficient of pi in the integral's result.
    The sum is S = sum_{j=0 to 25} C(50, 2j) * 4^j * (49!! * (2j-1)!!) / (50+2j)!!
    """
    total_sum = Fraction(0)
    d_fact_49 = double_factorial(49)

    for j in range(26):
        k = 2 * j
        
        # Binomial coefficient C(50, 2j)
        comb = combinations(50, k)
        
        # Power of 4
        pow4 = 4**j
        
        # Double factorials
        d_fact_k_minus_1 = double_factorial(k - 1)
        d_fact_50_plus_k = double_factorial(50 + k)
        
        # Numerator and denominator for this term
        numerator = comb * pow4 * d_fact_49 * d_fact_k_minus_1
        denominator = d_fact_50_plus_k
        
        term = Fraction(numerator, denominator)
        total_sum += term
        
    print("The final sum is:")
    # The final equation is I = pi * total_sum
    # Let's print the parts for the j=0 and j=1 terms as an example.
    # j=0:
    j=0
    k=0
    comb = combinations(50, k)
    pow4 = 4**j
    d_fact_k_minus_1 = double_factorial(k - 1)
    d_fact_50_plus_k = double_factorial(50 + k)
    # The term is C(50,0)*4^0 * (49!!*(-1)!!)/(50)!! = 1 * 1 * 49!! / 50!!
    print(f"Term for j=0 is C(50, 0) * 4^0 * (49!! * (-1)!!) / (50)!! which is {Fraction(comb * pow4 * d_fact_49 * d_fact_k_minus_1, d_fact_50_plus_k)}")

    # j=1:
    j=1
    k=2
    comb = combinations(50, k)
    pow4 = 4**j
    d_fact_k_minus_1 = double_factorial(k - 1)
    d_fact_50_plus_k = double_factorial(50 + k)
    # The term is C(50,2)*4^1 * (49!!*1!!)/(52)!!
    print(f"Term for j=1 is C(50, 2) * 4^1 * (49!! * 1!!) / (52)!! which is {Fraction(comb * pow4 * d_fact_49 * d_fact_k_minus_1, d_fact_50_plus_k)}")

    print(f"And so on... the total sum evaluates to the fraction:")
    print(total_sum)

main()