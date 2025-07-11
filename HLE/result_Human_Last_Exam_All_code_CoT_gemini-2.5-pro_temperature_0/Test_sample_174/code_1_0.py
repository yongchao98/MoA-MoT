from fractions import Fraction

# The 6th Bernoulli number, B_6, is 1/42.
B6 = Fraction(1, 42)

# Genus of the curve
g = 3

# Step 1: Calculate the orbifold Euler characteristic of M_3
# Formula: chi(M_g) = B_{2g} / (2g * (2g - 2))
# For g=3, chi(M_3) = B_6 / (6 * 4) = B_6 / 24
chi_M3_num = B6.numerator
chi_M3_den = B6.denominator * 24
chi_M3 = Fraction(chi_M3_num, chi_M3_den)

# Step 2: Calculate the orbifold Euler characteristic of H_3
# Formula: chi(H_g) = (1/2) * zeta(1-2g) = (1/2) * (-B_{2g} / (2g))
# For g=3, chi(H_3) = -B_6 / (2 * 6) = -B_6 / 12
chi_H3_num = -B6.numerator
chi_H3_den = B6.denominator * 12
chi_H3 = Fraction(chi_H3_num, chi_H3_den)

# Step 3: Calculate the final result by subtraction
# chi([U/G]) = chi(M_3) - chi(H_3)
result = chi_M3 - chi_H3

# Print the calculation step-by-step
print(f"The genus of a smooth plane quartic is g = (4-1)*(4-2)/2 = {g}.")
print(f"The 6th Bernoulli number is B_6 = {B6}.")
print(f"The orbifold Euler characteristic of the moduli stack of genus 3 curves is:")
print(f"chi(M_3) = B_6 / (2*g * (2*g - 2)) = {B6} / (6 * 4) = {B6} / 24 = {chi_M3}")
print(f"The orbifold Euler characteristic of the hyperelliptic locus of genus 3 is:")
print(f"chi(H_3) = -B_6 / (2 * g) = -({B6}) / (2 * 3) = -({B6}) / 12 = {chi_H3}")
print(f"The orbifold Euler characteristic of the quotient stack [U/G] is the difference:")
print(f"chi([U/G]) = chi(M_3) - chi(H_3) = {chi_M3} - ({chi_H3}) = {result}")
print("\nFinal Answer:")
print(f"{chi_M3} - ({chi_H3}) = {result}")
