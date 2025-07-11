import math
from fractions import Fraction

def calculate_residue(n):
    """
    Calculates the residue of f(z) at z = -n for n>=1.
    Res(f, -n) = (2n * (-1)^n) / ((2n+3) * n!)
    """
    if n == 0:
        return Fraction(0)
    
    numerator = 2 * n * ((-1)**n)
    denominator = (2 * n + 3) * math.factorial(n)
    return Fraction(numerator, denominator)

# Calculate residues for poles inside the contours
res_n1 = calculate_residue(1)
res_n2 = calculate_residue(2)
res_n3 = calculate_residue(3)

# Winding numbers for C1
ind_c1 = {-1: 1, -2: -1, -3: 1}

# Winding numbers for C2
ind_c2 = {-1: 1, -2: 1, -3: 0}

# Total contribution from each pole (Sum of (Ind_C1(zk) + Ind_C2(zk)) * Res(f, zk))
poles_to_consider = [-1, -2, -3]
residues = {-1: res_n1, -2: res_n2, -3: res_n3}
total_sum_of_residues = Fraction(0)

print("Sum of Integrals = 2 * pi * i * [ Sum of (Ind_C1(zk) + Ind_C2(zk)) * Res(f, zk) ]")
print("="*60)
full_calculation_str = "2 * pi * i * ( "

# Pole at z = -1
total_ind_n1 = ind_c1[-1] + ind_c2[-1]
term1_str = f"({ind_c1[-1]} + {ind_c2[-1]}) * ({res_n1})"
full_calculation_str += term1_str
total_sum_of_residues += total_ind_n1 * res_n1

# Pole at z = -2
total_ind_n2 = ind_c1[-2] + ind_c2[-2]
term2_str = f" + ({ind_c1[-2]} + {ind_c2[-2]}) * ({res_n2})"
full_calculation_str += term2_str
total_sum_of_residues += total_ind_n2 * res_n2

# Pole at z = -3
total_ind_n3 = ind_c1[-3] + ind_c2[-3]
term3_str = f" + ({ind_c1[-3]} + {ind_c2[-3]}) * ({res_n3})"
full_calculation_str += term3_str
full_calculation_str += " )"
total_sum_of_residues += total_ind_n3 * res_n3

print(f"Calculation breakdown: \n{full_calculation_str}\n")

sum_expr_str = f"2 * pi * i * ( {total_ind_n1} * ({res_n1}) + {total_ind_n2} * ({res_n2}) + {total_ind_n3} * ({res_n3}) )"
print(f"Simplified sum of contributions: \n{sum_expr_str}\n")

term1_val = total_ind_n1 * res_n1
term2_val = total_ind_n2 * res_n2
term3_val = total_ind_n3 * res_n3
final_sum_str = f"2 * pi * i * ( {term1_val} + {term2_val} + {term3_val} )"
print(f"Value of contributions: \n{final_sum_str}\n")

total_sum_fraction = term1_val + term2_val + term3_val
final_total_str = f"2 * pi * i * ( {total_sum_fraction} )"
print(f"Total sum of residues term: \n{final_total_str}\n")

final_integral_val_coeff = 2 * total_sum_fraction
final_integral_str = f"i * ( {final_integral_val_coeff.numerator}*pi / {final_integral_val_coeff.denominator} )"
print(f"The sum of the integrals is: \n{final_integral_str}\n")

imaginary_part_str = f"{final_integral_val_coeff.numerator}*pi / {final_integral_val_coeff.denominator}"
print(f"The imaginary part of the sum of the integrals is: {imaginary_part_str}")

final_numerical_value = float(final_integral_val_coeff) * math.pi
print(f"Numerical value: {final_numerical_value}")