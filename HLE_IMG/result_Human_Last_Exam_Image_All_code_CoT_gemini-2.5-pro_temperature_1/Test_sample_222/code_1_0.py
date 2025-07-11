import math
from fractions import Fraction

def calculate_residue(n):
    """Calculates the residue of f(z) at the pole z = -n."""
    # Using floating point numbers for the final numerical result
    numerator = 2 * n
    denominator = 2 * n + 3
    sign = (-1)**n
    factorial_n = math.factorial(n)
    residue = (numerator / denominator) * (sign / factorial_n)
    return residue

# Calculate the required residues for poles z=-2 and z=-3
res_f_minus_2 = calculate_residue(2)
res_f_minus_3 = calculate_residue(3)

# The imaginary part of the sum of integrals is 2 * pi * (2 * Res(f, -2) - Res(f, -3))
term_in_parenthesis = 2 * res_f_minus_2 - res_f_minus_3
imaginary_part = 2 * math.pi * term_in_parenthesis

# Using Fractions for exact representation in the explanation
res_f_minus_2_frac = Fraction(2 * 2, 2 * 2 + 3) * Fraction((-1)**2, math.factorial(2))
res_f_minus_3_frac = Fraction(2 * 3, 2 * 3 + 3) * Fraction((-1)**3, math.factorial(3))
term_in_parenthesis_frac = 2 * res_f_minus_2_frac - res_f_minus_3_frac
final_expression_frac = 2 * term_in_parenthesis_frac

print("The sum of integrals, I, simplifies to I = 2 * pi * i * [2 * Res(f, -2) - Res(f, -3)]")
print("\nCalculating the residues:")
print(f"Res(f, -2) = (2*({2}))/(2*({2})+3) * ((-1)^({2}))/({2}!) = {res_f_minus_2_frac.numerator}/{res_f_minus_2_frac.denominator}")
print(f"Res(f, -3) = (2*({3}))/(2*({3})+3) * ((-1)^({3}))/({3}!) = {res_f_minus_3_frac.numerator}/{res_f_minus_3_frac.denominator}")
print("\nSubstituting the residues into the expression for I:")
print(f"I = 2 * pi * i * [2 * ({res_f_minus_2_frac.numerator}/{res_f_minus_2_frac.denominator}) - ({res_f_minus_3_frac.numerator}/{res_f_minus_3_frac.denominator})]")
print(f"I = 2 * pi * i * [{2 * res_f_minus_2_frac} - ({res_f_minus_3_frac})]")
print(f"I = 2 * pi * i * [{term_in_parenthesis_frac.numerator}/{term_in_parenthesis_frac.denominator}]")
print("\nThe imaginary part of the sum is therefore:")
print(f"Im(I) = 2 * pi * ({term_in_parenthesis_frac.numerator}/{term_in_parenthesis_frac.denominator}) = ({final_expression_frac.numerator}*pi)/{final_expression_frac.denominator}")
print("\nThe final numerical value is:")
print(imaginary_part)