import math
from fractions import Fraction

# Step 1: Define the function and its poles.
# f(z) = z * Gamma(z) / (z - 3/2) = Gamma(z+1) / (z - 3/2)
# Poles are at z=3/2 and z = -1, -2, -3, ...

# Step 2: Calculate residues at the relevant poles.
# Residue at a simple pole z_0 is lim_{z->z_0} (z-z_0)f(z).
# For poles from Gamma function at z = -n (n=1,2,3...), the residue of Gamma(z+1) is (-1)^(n-1)/(n-1)!
# Res(f, z_0) = Res(Gamma(z+1), z_0) / (z_0 - 3/2)

# Res(f, 3/2) = Gamma(3/2 + 1) = Gamma(5/2) = (3/2)*(1/2)*sqrt(pi) = 3*sqrt(pi)/4
# This has a sqrt(pi) term. We'll track it separately.
res_1_5_factor = Fraction(3, 4)

# Res(f, -1): n=1. Res = ((-1)^0 / 0!) / (-1 - 3/2) = 1 / (-5/2) = -2/5
res_neg_1 = Fraction(-2, 5)

# Res(f, -2): n=2. Res = ((-1)^1 / 1!) / (-2 - 3/2) = -1 / (-7/2) = 2/7
res_neg_2 = Fraction(2, 7)

# Res(f, -3): n=3. Res = ((-1)^2 / 2!) / (-3 - 3/2) = (1/2) / (-9/2) = -1/9
res_neg_3 = Fraction(-1, 9)

# Step 3 & 4: Determine winding numbers and apply the residue theorem.
# C2 (orange): CCW, encloses 1.5, -1, -2. Winding numbers are +1.
# C1 (blue): a series of loops. From the arrows, the loop directions alternate.
# The loop containing (0,1) is CCW.
# Thus, loop containing (1,2) [and pole 1.5] is CW (W=-1).
# Loop containing (-1,0) [and pole -1] is CW (W=-1).
# Loop containing (-2,-1) [and pole -2] is CCW (W=+1).
# Loop containing (-3,-2) [and pole -3] is CW (W=-1).
# The drawing for C1 seems to stop before enclosing -4.

# Step 5: Sum the integrals by summing the winding numbers for each pole.
# Total Integral I = 2*pi*i * sum over poles k of [ (W(C1,k) + W(C2,k)) * Res(f,k) ]

# For pole z = 1.5:
w_sum_1_5 = -1 + 1  # W(C1, 1.5) + W(C2, 1.5)
# For pole z = -1:
w_sum_neg_1 = -1 + 1  # W(C1, -1) + W(C2, -1)
# For pole z = -2:
w_sum_neg_2 = 1 + 1   # W(C1, -2) + W(C2, -2)
# For pole z = -3:
w_sum_neg_3 = -1 + 0  # W(C1, -3) + W(C2, -3)

# The sum contains contributions only from poles where the sum of winding numbers is non-zero.
# The sqrt(pi) term from Res(f, 1.5) cancels out, as does the term from Res(f, -1).

# Step 6: Calculate the final sum.
total_residue_contribution = (w_sum_neg_2 * res_neg_2) + (w_sum_neg_3 * res_neg_3)

# Step 7: Print the calculation and the imaginary part.
print("The sum of the integrals is I = 2*pi*i * sum[ (W(C1,k) + W(C2,k)) * Res(f,k) ]")
print("We only need to consider poles where the sum of winding numbers is not zero.")
print("These are the poles at z = -2 and z = -3.\n")

print("For z = -2:")
print(f"  Res(f, -2) = {res_neg_2.numerator}/{res_neg_2.denominator}")
print(f"  Sum of winding numbers W_sum = W(C1,-2) + W(C2,-2) = 1 + 1 = {w_sum_neg_2}")
print(f"  Contribution = {w_sum_neg_2} * ({res_neg_2.numerator}/{res_neg_2.denominator}) = {w_sum_neg_2 * res_neg_2.numerator}/{res_neg_2.denominator}\n")

print("For z = -3:")
print(f"  Res(f, -3) = {res_neg_3.numerator}/{res_neg_3.denominator}")
print(f"  Sum of winding numbers W_sum = W(C1,-3) + W(C2,-3) = -1 + 0 = {w_sum_neg_3}")
print(f"  Contribution = {w_sum_neg_3} * ({res_neg_3.numerator}/{res_neg_3.denominator}) = {w_sum_neg_3 * res_neg_3.numerator}/{res_neg_3.denominator}\n")

print(f"Total contribution from residues = (2 * 2/7) + (-1 * -1/9) = 4/7 + 1/9")
print(f"4/7 + 1/9 = {Fraction(4,7) + Fraction(1,9).numerator}/{Fraction(4,7).denominator * Fraction(1,9).denominator} = {total_residue_contribution.numerator}/{total_residue_contribution.denominator}")

print(f"\nSo, the sum of integrals is I = 2*pi*i * ({total_residue_contribution.numerator}/{total_residue_contribution.denominator})")
final_coefficient_numerator = 2 * total_residue_contribution.numerator
final_coefficient_denominator = total_residue_contribution.denominator
print(f"I = i * ({final_coefficient_numerator}*pi / {final_coefficient_denominator})")

imaginary_part = (final_coefficient_numerator * math.pi) / final_coefficient_denominator
print(f"\nThe imaginary part of the sum is {final_coefficient_numerator}*pi/{final_coefficient_denominator} which is approximately {imaginary_part}")

final_answer = (86 * math.pi) / 63
# The prompt asks for the final answer in a specific format.
# print(f"<<<{final_answer}>>>")