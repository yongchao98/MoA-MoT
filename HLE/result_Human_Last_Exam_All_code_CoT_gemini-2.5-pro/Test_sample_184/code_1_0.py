import math
from fractions import Fraction
import sympy

def solve():
    """
    Solves the problem by calculating the sum of 1/n^2 for n in the specified set S.
    """
    
    # Step 1: Explain the problem setup using the Euler product.
    print("Let S be the set of positive integers n such that for any prime p, the exponent of p in the prime factorization of n is not 1, 2, or 5.")
    print("We want to evaluate Sum_{n in S} (1/n^2).")
    print("Since the condition on n is multiplicative, the sum can be written as an Euler product over all primes p:")
    print("Sum = Product_p [ Sum_{e in E} (1 / (p^e)^2) ]")
    print("where E is the set of allowed exponents E = {0, 3, 4, 6, 7, 8, ...} = N_0 \\ {1, 2, 5}.\n")

    # Step 2: Connect to the Frobenius Coin Problem.
    print("Step 2: Connect the set of exponents to the Frobenius Coin Problem.")
    print("The set of numbers that can be written in the form 3a + 4b for non-negative integers a, b is N_0 \\ {1, 2, 5}.")
    print("This means the set of allowed exponents E is precisely the set of values {3a + 4b | a >= 0, b >= 0}.")
    print("This allows us to rewrite the inner sum over exponents.\n")

    # Step 3: Simplify the per-prime sum.
    print("Step 3: Simplify the sum for each prime p.")
    print("Let x = p^(-2). The sum for a prime p is Sum_{e in E} x^e.")
    print("Using the Frobenius coin result, this sum becomes:")
    print("Sum_{a>=0, b>=0} x^(3a+4b) = (Sum_{a>=0} x^(3a)) * (Sum_{b>=0} x^(4b))")
    print("These are two geometric series:")
    print("= [1 / (1 - x^3)] * [1 / (1 - x^4)]")
    print("Substituting x = p^(-2) back, the term for each prime p is:")
    print("1 / ( (1 - p^-6) * (1 - p^-8) )\n")

    # Step 4: Express the total sum in terms of the Riemann zeta function.
    print("Step 4: Express the total sum using the Riemann zeta function.")
    print("The total sum is the product of these terms over all primes:")
    print("Sum = Product_p [ 1 / ( (1 - p^-6) * (1 - p^-8) ) ]")
    print("= (Product_p [1 / (1 - p^-6)]) * (Product_p [1 / (1 - p^-8)])")
    print("This is the product of two Riemann zeta functions, based on the Euler product formula zeta(s) = Product_p [1 / (1 - p^-s)].")
    print("Sum = zeta(6) * zeta(8)\n")

    # Step 5: Calculate the final value.
    print("Step 5: Calculate the final value using known formulas for zeta(2k).")
    # zeta(6) = pi^6 / 945
    # zeta(8) = pi^8 / 9450
    zeta6_num = 1
    zeta6_den = 945
    
    zeta8_num = 1
    zeta8_den = 9450

    # Represent as fractions for precise arithmetic
    frac6 = Fraction(zeta6_num, zeta6_den)
    frac8 = Fraction(zeta8_num, zeta8_den)
    
    final_frac = frac6 * frac8
    
    # Using sympy to display pi nicely
    pi = sympy.symbols('pi')
    
    print(f"We know that zeta(6) = pi^6 / {zeta6_den} and zeta(8) = pi^8 / {zeta8_den}.")
    print(f"So, the sum is (pi^6 / {zeta6_den}) * (pi^8 / {zeta8_den}).")
    print(f"Sum = ({final_frac.numerator} * pi^14) / {final_frac.denominator}")
    
    final_answer_str = f"({final_frac.numerator}*pi^14)/{final_frac.denominator}"
    return final_answer_str

final_answer = solve()
print(f"\nFinal Answer: {final_answer}")

# Format the final answer as requested
print(f'<<<{final_answer}>>>')
