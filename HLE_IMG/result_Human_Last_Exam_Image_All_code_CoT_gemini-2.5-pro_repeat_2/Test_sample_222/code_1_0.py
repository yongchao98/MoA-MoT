import math
from fractions import Fraction

def calculate_residue(n):
    """
    Calculates the residue of f(z) at the pole z = -n.
    Res(f, -n) = (n / (n + 3/2)) * ((-1)^n / n!)
    """
    numerator = n
    denominator = Fraction(n) + Fraction(3, 2)
    term1 = Fraction(numerator) / denominator
    
    term2_numerator = (-1)**n
    term2_denominator = math.factorial(n)
    term2 = Fraction(term2_numerator, term2_denominator)
    
    return term1 * term2

# From the analysis, the sum of integrals is I = 2*pi*i * (2*Res(f, -1) + Res(f, -3))
# Calculate the required residues
res_minus_1 = calculate_residue(1)
res_minus_3 = calculate_residue(3)

# The sum of residues term S = 2*Res(f,-1) + Res(f,-3)
S = 2 * res_minus_1 + res_minus_3

# The total integral is I = 2 * pi * i * S
# The imaginary part is Im(I) = 2 * pi * S
imaginary_part_coeff = 2 * S

print("The sum of the integrals is given by I = 2*pi*i * (2 * Res(f, -1) + Res(f, -3)).")
print(f"Calculating the residues:")
print(f"Res(f, -1) = {res_minus_1.numerator}/{res_minus_1.denominator}")
print(f"Res(f, -3) = {res_minus_3.numerator}/{res_minus_3.denominator}")
print("\nSubstituting the values into the formula for the imaginary part:")
print(f"Im(I) = 2 * pi * (2 * ({res_minus_1.numerator}/{res_minus_1.denominator}) + ({res_minus_3.numerator}/{res_minus_3.denominator}))")
print(f"Im(I) = 2 * pi * (({2 * res_minus_1.numerator}/{res_minus_1.denominator}) - {abs(res_minus_3.numerator)}/{res_minus_3.denominator})")
print(f"Im(I) = 2 * pi * ({S.numerator}/{S.denominator})")
print(f"Im(I) = ({imaginary_part_coeff.numerator} * pi) / {imaginary_part_coeff.denominator}")

# Calculate the final numerical value
numerical_value = float(imaginary_part_coeff) * math.pi
print(f"\nThe numerical value of the imaginary part is approximately: {numerical_value}")